"""
    (x::SupplyChainEnv)(action::Vector{T} where T <: Real)

Apply an `action` (replenishment requests) on the `SupplyChainEnv` and step
forward one simulation period.
"""
function (x::SupplyChainEnv)(action::Vector{T} where T <: Real)
    #validate action input
    @assert all(action .>= 0) "Reorder actions cannot be negative."
    @assert length(action) == length(x.materials)*ne(x.network) "Reorder action vector must have length num_products * num_edges."
    
    #increase period counter
    x.period += 1

    #find previous period's inventories
    prev_inv_on_hand = filter(:period => j -> j == x.period-1, x.inv_on_hand, view=true) #previous inventory level
    prev_onhand_grp = groupby(prev_inv_on_hand, [:node, :material])
    prev_inv_pipeline = filter(:period => j -> j == x.period-1, x.inv_pipeline, view=true) #previous inventory level
    prev_pipeln_grp = groupby(prev_inv_pipeline, [:arc, :material])
    
    #intialize next period on-hand and pipeline inventories with previous inventories
    for p in x.materials
        for n in vertices(x.network)
            prev = prev_onhand_grp[(node = n, material = p)][1,:level] #previous on hand inventory
            push!(x.inv_on_hand, [x.period, n, p, prev, 0]) #intialize with previous inventory levels
        end
        for a in edges(x.network)
            prev = prev_pipeln_grp[(arc = (a.src,a.dst), material = p)][1,:level] #previous pipeline inventory
            push!(x.inv_pipeline, [x.period, (a.src,a.dst), p, prev]) #intialize with previous inventory levels
        end
    end
    
    #move active shipments forward one period
    x.shipments.lead .-= 1

    #move active production forward one period
    x.production.lead .-= 1

    #reshape action
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    act = reshape(action, (length(x.materials), length(arcs)))

    #place requests
    place_requests!(x, act, arcs)

    #update production and associated pipeline inventories
    update_production!(x)

    #update on hand and pipeline inventories due to arrived shipments
    arrivals = update_shipments!(x)

    #discard any excess inventory
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)

    #markets open and demand occurs
    simulate_markets!(x)

    #update inventory positions at each node
    update_positions!(x)

    if x.options[:evaluate_profit]
        #calculate profit at each node
        calculate_profit!(x, arrivals)
        #update reward (current profit). Should be revised for RL
        x.reward = sum(filter(:period => j -> j == x.period, x.profit, view=true).value)
    end
end

"""
    place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)

Place inventory replenishment requests throughout the network.
"""
function place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)

    #exit if no action
    if iszero(act)
        for a in arcs, p in x.materials
            push!(x.replenishments, [x.period, a, p, 0, 0, 0, 0, missing])
        end
        return
    end

    #extract info
    mats = x.materials
    bom = x.bill_of_materials

    #store copy of action vector
    requests = copy(act)

    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = deepcopy(Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers))

    #sample lead times
    leads = Dict((a,p) => rand(get_prop(x.network, a, :lead_time)[p]) for a in edges(x.network), p in mats)

    #non source nodes
    nonsources = [n for n in vertices(x.network) if !isempty(inneighbors(x.network, n))]

    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply

    #get backlog
    if x.options[:backlog] && x.period > 1
        last_orders_df = filter(:period => i -> i == x.period-1, x.replenishments, view=true) #replenishment orders from previous period
        last_orders_grp = groupby(last_orders_df, :material) #group by material
    end

    #place requests
    for (i,p) in enumerate(mats) #loop by materials
        for n in nonsources #loop by nodes placing requests
            sup_priority = get_prop(x.network, n, :supplier_priority)[p] #get supplier priority list
            for (l, sup) in enumerate(sup_priority) #loop by supplier priority
                a = (sup, n) #arc
                j = findfirst(k -> k == a, arcs) #find index for that arc in the action matrix
                
                #amount requested
                amount = act[i,j] 
                if x.options[:backlog] && x.period > 1 #add previous period's backlog
                    if x.options[:reallocate] && l == 1
                        amount += filter(:arc => k -> k == (sup_priority[end],n), last_orders_grp[(material = p,)], view=true).unfulfilled[1]
                    elseif !x.options[:reallocate]
                        amount += filter(:arc => k -> k == a, last_orders_grp[(material = p,)], view=true).unfulfilled[1]
                    end
                end

                #continue to next iteration if request is 0
                if iszero(amount) 
                    push!(x.replenishments, [x.period, a, p, requests[i,j], 0, 0, 0, missing])
                    continue 
                end

                #available supply
                supply = supply_grp[(node = a[1], material = p)].level[1]

                #accept or adjust requests
                if a[1] in x.producers
                    #try to satisfy with on hand inventory first
                    accepted_inv = min(amount, supply) 
                    amount -= accepted_inv #get new amount pending

                    #try to satisfy with by-products
                    sched_supply = filter([:arc, :material] => (k1, k2) -> k1 == (a[1],a[1]) && k2 == p, x.production)[:,[:amount,:lead]] #check if any inventory is scheduled for production, but not commited to downstream node
                    sort!(sched_supply, :lead) #sort to use scheduled supply that is closest to finishing
                    sched_supp, sched_lead = isempty(sched_supply) ? [0,0] : [sched_supply.amount, sched_supply.lead]
                    accepted_sched_k = [] #try to satisfy with scheduled production that is not commited
                    for k in sched_supp
                        accepted_k = min(amount, k) #satisfy if possible with scheduled supply
                        push!(accepted_sched_k, accepted_k)
                        amount -= accepted_k #update pending amount
                    end
                    accepted_sched = sum(accepted_sched_k)

                    #try to satisfy with material production
                    capacity = [] #get production capacity
                    mat_supply = [] #store max capacity based on raw material consumption for each raw material
                    imats = findall(k -> !iszero(k), bom[:,i]) #indices of materials involved with production of p
                    if !isempty(imats) #only add capacity if there is a non-zero bom for that material
                        push!(capacity, capacities[a[1]][p])
                    else
                        push!(capacity, 0)
                    end
                    for ii in imats
                        sup_pp = supply_grp[(node = a[1], material = mats[ii])].level[1] #supply of material involved in BOM
                        if bom[ii,i] < 0 #only account for raw materials that are in the BOM
                            push!(mat_supply, - sup_pp / bom[ii,i])
                        elseif bom[ii,i] > 0 #add capacity constraint for any co-products (scaled by stoichiometry)
                            push!(capacity, capacities[a[1]][mats[ii]] / bom[ii,i])
                        end
                    end
                    accepted_prod = min(amount, capacity..., mat_supply...) #fulfill remaining with available capacity
                    amount -= accepted_prod #update pending amount
                    capacities[a[1]][p] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)

                    #total accepted request
                    accepted = accepted_inv + accepted_sched + accepted_prod

                    #schedule production
                    prod_time = get_prop(x.network, a[1], :production_time)[p] #production time
                    accepted_prod > 0 && push!(x.production, (a, p, accepted_prod, prod_time))

                    #update commited scheduled production of byproducts to new destination (decrease amount sent to self storage and add commited amount to arc)
                    for (k, acpt_sched) in enumerate(accepted_sched_k)
                        if acpt_sched > 0
                            x.production[(string.(x.production.arc) .== string((a[1],a[1]))) .&
                                         (x.production.material .== p) .&
                                         (x.production.lead .== sched_lead[k]) .&
                                         (x.production.amount .== sched_supp[k]), :amount] .-= acpt_sched
                            push!(x.production, (a, p, acpt_sched, sched_lead[k]))
                        end
                    end

                    #update inventories
                    if accepted_inv > 0
                        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                      (x.inv_on_hand.node .== a[1]) .&
                                      (x.inv_on_hand.material .== p), :level] .-= accepted_inv
                        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                                       (x.inv_pipeline.material .== p), :level] .+= accepted_inv
                    end
                    if accepted_prod > 0
                        for ii in imats
                            if bom[ii,i] < 0 #reactant consumed
                                x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                              (x.inv_on_hand.node .== a[1]) .&
                                              (x.inv_on_hand.material .== mats[ii]), :level] .+= accepted_prod * bom[ii,i]
                            elseif bom[ii,i] > 0 #coproducts scheduled for production
                                a1 = (a[1],a[1]) #production of byproduct is going to inventory holding at current node (not being shipped)
                                push!(x.production, (a1, mats[ii], accepted_prod*bom[ii,i], prod_time))
                                capacities[a[1]][mats[ii]] -= accepted_prod*bom[ii,i]
                            end
                        end
                    end
                else
                    accepted = min(amount, supply) #accepted request
                    if accepted > 0
                        amount -= accepted #update remaining amount
                        #update inventories (on hand and pipeline)
                        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                    (x.inv_on_hand.node .== a[1]) .&
                                    (x.inv_on_hand.material .== p), :level] .-= accepted
                        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                                    (string.(x.inv_pipeline.arc) .== string(a)) .&
                                    (x.inv_pipeline.material .== p), :level] .+= accepted
                    end
                end

                #check if some was not accepted and reallocate
                new_alloc = missing
                if amount > 0 #reallocate unfulfilled request to next priority supplier
                    if x.options[:reallocate]
                        next_sup = findfirst(k -> k == a[1], sup_priority) + 1 #get next in line
                        if next_sup <= length(sup_priority) #check that there is a next one in line
                            new_sup = sup_priority[next_sup] #next supplier in line
                            jj = findfirst(k -> k == (new_sup, a[2]), arcs) #find index for that arc in the action matrix
                            act[i,jj] += amount #add unfulfilled to the next supplier in the line
                            new_alloc = (new_sup, a[2]) #store new arc where reallocated
                        end
                    end
                end

                #store shipments and lead times
                if accepted > 0 #update active shipments
                    lead = leads[Edge(a[1],a[end]), p]
                    if a[1] in x.producers
                        if accepted_inv > 0
                            push!(x.shipments, [a, p, accepted_inv, lead])
                            push!(x.replenishments, [x.period, a, p, requests[i,j], accepted_inv, lead, amount, new_alloc])
                        end
                        if accepted_sched > 0
                            for (k, acpt_sched) in enumerate(accepted_sched_k)
                                push!(x.shipments, [a, p, acpt_sched, lead + sched_lead[k]])
                                push!(x.replenishments, [x.period, a, p, requests[i,j], acpt_sched, lead + sched_lead[k], amount, new_alloc])
                            end
                        end
                        if accepted_prod > 0
                            push!(x.shipments, [a, p, accepted_prod, lead + prod_time])
                            push!(x.replenishments, [x.period, a, p, requests[i,j], accepted_prod, lead + prod_time, amount, new_alloc])
                        end
                    else
                        push!(x.shipments, [a, p, accepted, lead])
                        push!(x.replenishments, [x.period, a, p, requests[i,j], accepted, lead, amount, new_alloc])
                    end
                else #no request made
                    push!(x.replenishments, [x.period, a, p, requests[i,j], 0, 0, amount, new_alloc])
                end
            end
        end
    end

    #updated production capacities (after commited production)
    for n in x.producers
        set_prop!(x.network, n, :production_capacity, capacities[n])
    end
end

"""
    update_production!(x::SupplyChainEnv)

Update completed production at each producer node and send to the arc that it
    was commited to.
"""
function update_production!(x::SupplyChainEnv)
    #extract info
    mats = x.materials
    bom = x.bill_of_materials

    #find production that has completed
    produced = filter(:lead => i -> i <= 0, x.production) #find active production that has completed
    for i in 1:nrow(produced) #add produced material to pipeline
        a, p, amount = produced[i, [:arc, :material, :amount]]
        #restore production capacities
        capacity = get_prop(x.network, a[1], :production_capacity)
        capacity[p] += amount
        j = findfirst(k -> k==p, mats)
        imats = findall(k -> k > 0, bom[:,j])
        for ii in imats
            capacity[mats[ii]] += amount*bom[ii,j]
        end
        set_prop!(x.network, a[1], :production_capacity, capacity)
        #update inventories
        if a[1] != a[2] #send produced material down pipeline
            x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                           (string.(x.inv_pipeline.arc) .== string(a)) .&
                           (x.inv_pipeline.material .== p), :level] .+= amount
        else #send produced material to storage at that node
            x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                          (x.inv_on_hand.node .== a[1]) .&
                          (x.inv_on_hand.material .== p), :level] .+= amount
        end
    end
    filter!(:lead => i -> i > 0, x.production) #remove produced inventory that was shipped (and any zero values)
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    arrivals = filter(:lead => i -> i <= 0, x.shipments) #find active shipments with 0 lead time
    for i in 1:nrow(arrivals)
        a, p, amount = arrivals[i, [:arc, :material, :amount]]
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== a[end]) .&
                      (x.inv_on_hand.material .== p), :level] .+= amount
        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                       (x.inv_pipeline.material .== p), :level] .-= amount
    end
    filter!(:lead => i -> i > 0, x.shipments) #remove shipments that arrived (and any zero values)
    return arrivals
end

"""
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventories
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material for easier lookup
    for n in vertices(x.network)
        node_max_inv = get_prop(x.network, n, :inventory_capacity)
        for p in x.materials
            max_inv = node_max_inv[p]
            isinf(max_inv) && continue
            onhand = onhand_grp[(node = n, material = p)].level[1]
            if onhand > max_inv
                x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                            (x.inv_on_hand.node .== n) .&
                            (x.inv_on_hand.material .== p), :level] .= max_inv
                x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                            (x.inv_on_hand.node .== n) .&
                            (x.inv_on_hand.material .== p), :discarded] .+= onhand - max_inv
            end
        end
    end
end


"""
    update_positions!(x::SupplyChainEnv)

Update inventory position and inventory level for all materials and nodes.
"""
function update_positions!(x::SupplyChainEnv)
    #filter data
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, :material) #group by material
    dmnd_df = filter(:period => i -> i == x.period, x.demand, view=true) #demand at markets
    dmnd_grp = groupby(dmnd_df, [:node, :material]) #group by node and material
    orders_df = filter([:period, :reallocated] => (i1, i2) -> i1 == x.period && ismissing(i2), x.replenishments, view=true) #replenishment orders from current period
    orders_grp = groupby(orders_df, :material) #group by materil

    for n in vertices(x.network), p in x.materials
        making = sum(filter([:arc, :material] => (j1, j2) -> j1[end] == n && j2 == p, x.production, view=true).amount) #commited production order
        upstream = sum(filter(:arc => j -> j[end] == n, pipeline_grp[(material = p,)], view=true).level) #in-transit inventory
        onhand = onhand_grp[(node = n, material = p)].level[1] #on_hand inventory
        backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
        backlog = 0 #initialize backlog for orders placed by successors
        if x.options[:backlog]
            if n in x.markets #find demand backlog
                backlog += dmnd_grp[(node = n, material = p)].unfulfilled[1]
            else #find any unfulfilled replenishment request that was not reallocated
                backlog += sum(filter(:arc => i -> i[1] == n, orders_grp[(material = p,)], view=true).unfulfilled)
            end
            backorder = sum(filter(:arc => i -> i[end] == n, orders_grp[(material = p,)], view=true).unfulfilled)
        end
        push!(x.inv_level, [x.period, n, p, onhand - backlog]) #update inventory
        push!(x.inv_position, [x.period, n, p, onhand + making + upstream + backorder - backlog]) #update inventory
    end
end

"""
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    last_dmnd = filter(:period => i -> i == x.period-1, x.demand, view=true)
    last_d_group = groupby(last_dmnd, [:node, :material])
    for n in x.markets, p in x.materials
        dmnd_seq = get_prop(x.network, n, :demand_sequence)[p]
        freq = Bernoulli(get_prop(x.network, n, :demand_frequency)[p]) #demand frequency
        dmnd = get_prop(x.network, n, :demand_distribution)[p] #demand distribution
        q = iszero(dmnd_seq) ? rand(freq) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
        demand = [x.period, n, p, q, 0., 0.] #initialize demand vector to store in df
        if x.options[:backlog] && x.period > 1 #add previous backlog to the quantity requested at the market
            q += last_d_group[(node = n, material = p)].unfulfilled[1]
        end
        if q > 0
            inv = onhand_grp[(node = n, material = p)].level[1]
            demand[5] = min(q, inv) #sales
            demand[6] = max(q - inv, 0.) #unfilfilled
        end
        push!(x.demand, demand) #update df

        #update end of period inventory (subtract sales)
        if demand[5] > 0
            x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                        (x.inv_on_hand.node .== n) .&
                        (x.inv_on_hand.material .== p), :level] .-= demand[5]
        end
    end
end

"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    on_hand_df = filter([:period, :level] => (j1, j2) -> j1 == x.period && !isinf(j2), x.inv_on_hand, view=true) #on_hand inventory 
    onhand_grp = groupby(on_hand_df, [:node, :material])
    orders_df = filter(:period => i -> i == x.period, x.replenishments, view=true) #replenishment orders
    orders_grp = groupby(orders_df, :material)
    sales_df = filter(:period => j -> j == x.period, x.demand, view=true) #sales at markets
    sales_grp = groupby(sales_df, [:node, :material])
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, [:arc, :material])
    for n in vertices(x.network)
        profit = 0 #initialize node profit
        for p in x.materials
            #get costs
            #holding cost
            hold_cost = get_prop(x.network, n, :holding_cost)[p]
            if hold_cost > 0
                onhand = onhand_grp[(node = n, material = p)].level
                if !isempty(onhand)
                    profit -= hold_cost * onhand[1]
                end
            end
            #production cost
            if n in x.producers
                prod_cost = get_prop(x.network, n, :production_cost)[p]
                if prod_cost > 0
                    produced = sum(filter([:arc] => j -> j[1] == n, orders_grp[(material = p,)], view=true).accepted)
                    profit -= prod_cost * produced
                end
            #sales profit at markets (and penalize for unfulfilled demand)
            elseif n in x.markets
                sales_price = get_prop(x.network, n, :sales_price)[p]
                dmnd_penalty = get_prop(x.network, n, :demand_penalty)[p]
                if sales_price > 0 || dmnd_penalty > 0
                    sold, unfilled = sales_grp[(node = n, material = p)][1, [:sold, :unfulfilled]]
                    profit += sales_price * sold
                    profit -= dmnd_penalty * unfilled
                end
            end
            #pay suppliers for received inventory and pay transportation cost
            for pred in inneighbors(x.network, n)
                price = get_prop(x.network, pred, n, :sales_price)[p]
                trans_cost = get_prop(x.network, pred, n, :transportation_cost)[p]
                if price > 0 #pay purchase of inventory
                    purchased = filter([:arc, :material] => (j1, j2) -> j1 == (pred, n) && j2 == p, arrivals, view=true).amount
                    if !isempty(purchased)
                        profit -= purchased[1] * price
                    end
                end
                if trans_cost > 0 #pay transportation (pipeline hoding) cost (assume it is paid to a third party)
                    intransit = pipeline_grp[(arc = (pred, n), material = p)].level[1]
                    profit -= intransit * trans_cost
                end
            end
            #receive payment for delivered inventory
            for succ in outneighbors(x.network, n)
                price = get_prop(x.network, n, succ, :sales_price)[p]
                if price > 0 #receive payment for delivered inventory
                    sold = filter([:arc, :material] => (j1, j2) -> j1 == (n, succ) && j2 == p, arrivals, view=true).amount
                    if !isempty(sold)
                        profit += sold[1] * price
                    end
                end
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end
