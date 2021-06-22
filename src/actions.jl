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

    #intialize next period on-hand and pipeline inventories with previous inventories
    for n in vertices(x.network), p in x.materials
        prev = filter(j -> j.period == x.period-1 && j.node == n && j.material == p, x.inv_on_hand).level[1] #previous inventory level
        push!(x.inv_on_hand, [x.period, n, p, prev, 0]) #intialize with previous inventory levels
    end
    for a in edges(x.network), p in x.materials
        prev = filter(j -> j.period == x.period-1 && j.arc == (a.src,a.dst) && j.material == p, x.inv_pipeline).level[1] #previous inventory level
        push!(x.inv_pipeline, [x.period, (a.src,a.dst), p, prev]) #intialize with previous inventory levels
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
    enforce_inventory_limits!(x)

    #update inventory positions for non-market nodes
    for n in union(x.producers, x.distributors), p in x.materials
        update_position!(x, n, p)
    end

    #markets open and demand occurs
    simulate_markets!(x)

    #calculate profit at each node
    calculate_profit!(x, arrivals)

    #update reward (current profit). Should be revised for RL
    x.reward = sum(filter(j -> j.period == x.period, x.profit).value)
end

"""
    place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)

Place inventory replenishment requests throughout the network.
"""
function place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)

    #extract info
    mats = x.materials
    bom = x.bill_of_materials

    #store copy of action vector
    requests = copy(act)

    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = deepcopy(Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers))

    #sample lead times
    leads = Dict((a,p) => rand(get_prop(x.network, a, :lead_time)[p]) for a in edges(x.network), p in mats)

    nonsources = [n for n in vertices(x.network) if !isempty(inneighbors(x.network, n))]
    for (i,p) in enumerate(mats) #loop by materials
        for n in nonsources #loop by nodes placing requests
            sup_priority = get_prop(x.network, n, :supplier_priority)[p] #get supplier priority list
            for (l, sup) in enumerate(sup_priority) #loop by supplier priority
                a = (sup, n) #arc
                j = findfirst(k -> k == a, arcs) #find index for that arc in the action matrix
                amount = act[i,j] #amount requested
                if x.backlog && x.period > 1 #add previous period's backlog
                    if x.reallocate && l == 1
                        amount += filter(k -> k.material == p && k.arc == (sup_priority[end],n) && k.period == x.period-1, x.replenishments).unfulfilled[1]
                    elseif !x.reallocate
                        amount += filter(k -> k.material == p && k.arc == a && k.period == x.period-1, x.replenishments).unfulfilled[1]
                    end
                end
                if iszero(amount) #continue to next iteration if request is 0
                    push!(x.replenishments, [x.period, a, p, requests[i,j], 0, 0, 0, 0])
                    continue 
                end
                supply = filter(k -> k.period == x.period && k.node == a[1] && k.material == p, x.inv_on_hand).level[1] #on_hand inventory at supplier
                sched_supply = filter(k -> k.arc == (a[1],a[1]) && k.material == p, x.production)[:,[:amount,:lead]] #check if any inventory is scheduled for production, but not commited to downstream node
                sort!(sched_supply, :lead) #sort to use scheduled supply that is closest to finishing
                sched_supp, sched_lead = isempty(sched_supply) ? [0,0] : [sched_supply.amount, sched_supply.lead]
                #accept or adjust requests
                if a[1] in x.producers
                    capacity = [] #get production capacity
                    mat_supply = [] #store max capacity based on raw material consumption for each raw material
                    imats = findall(k -> !iszero(k), bom[:,i]) #indices of materials involved with production of p
                    if !isempty(imats) #only add capacity if there is a non-zero bom for that material
                        push!(capacity, capacities[a[1]][p])
                    else
                        push!(capacity, 0)
                    end
                    for ii in imats
                        sup_pp = filter(k -> k.period == x.period && k.node == a[1] && k.material == mats[ii], x.inv_on_hand).level[1] #supply of pp
                        if bom[ii,i] < 0 #only account for raw materials that are in the BOM
                            push!(mat_supply, - sup_pp / bom[ii,i])
                        elseif bom[ii,i] > 0 #add capacity constraint for any co-products (scaled by stoichiometry)
                            push!(capacity, capacities[a[1]][mats[ii]] / bom[ii,i])
                        end
                    end
                    accepted_inv = min(amount, supply) #try to satisfy with on hand inventory first
                    amount -= accepted_inv #get new amount pending
                    accepted_sched_k = [] #try to satisfy with scheduled production that is not commited
                    for k in sched_supp
                        accepted_k = min(amount, k) #satisfy if possible with scheduled supply
                        push!(accepted_sched_k, accepted_k)
                        amount -= accepted_k #update pending amount
                    end
                    accepted_sched = sum(accepted_sched_k)
                    accepted_prod = min(amount, capacity..., mat_supply...) #fulfill remaining with available capacity
                    amount -= accepted_prod #update pending amount
                    capacities[a[1]][p] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
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
                    amount -= accepted #update remaining amount
                    #update inventories (on hand and pipeline)
                    x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                  (x.inv_on_hand.node .== a[1]) .&
                                  (x.inv_on_hand.material .== p), :level] .-= accepted
                    x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                                   (string.(x.inv_pipeline.arc) .== string(a)) .&
                                   (x.inv_pipeline.material .== p), :level] .+= accepted
                end

                #chekc if some was not accepted and reallocate
                # unfulfilled = amount - accepted #amount not accepted
                new_alloc = missing
                if amount > 0 #reallocate unfulfilled request to next priority supplier
                    # @warn "Replenishment request made in period $(x.period) by node $(a[2]) to node $(a[1]) for material $p was reduced from $amount to $accepted due to insufficient production capacity or insufficient inventory."
                    if x.reallocate
                        next_sup = findfirst(k -> k == a[1], sup_priority) + 1 #get next in line
                        if next_sup <= length(sup_priority) #check that there is a next one in line
                            new_sup = sup_priority[next_sup]
                            # @warn "Reallocating non-accepted request for node $(a[2]) for material $p to node $new_sup (amount = $unfulfilled)"
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
    produced = filter(i -> i.lead <= 0, x.production) #find active production that has completed
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
    filter!(i -> i.lead > 0, x.production) #remove produced inventory that was shipped (and any zero values)
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    arrivals = filter(i -> i.lead <= 0, x.shipments) #find active shipments with 0 lead time
    for i in 1:nrow(arrivals)
        a, p, amount = arrivals[i, [:arc, :material, :amount]]
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== a[end]) .&
                      (x.inv_on_hand.material .== p), :level] .+= amount
        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                       (x.inv_pipeline.material .== p), :level] .-= amount
    end
    filter!(i -> i.lead > 0, x.shipments) #remove shipments that arrived (and any zero values)
    return arrivals
end

"""
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    for n in vertices(x.network), p in x.materials
        max_inv = get_prop(x.network, n, :inventory_capacity)[p]
        onhand = filter(j -> j.period == x.period && j.node == n && j.material == p, x.inv_on_hand).level[1] #on_hand inventory
        if onhand > max_inv
            # @warn "$(onhand - max_inv) units of material $p discarded at node $n because maximum inventory was exceeded in period $(x.period)."
            x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                          (x.inv_on_hand.node .== n) .&
                          (x.inv_on_hand.material .== p), :level] .= max_inv
            x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                          (x.inv_on_hand.node .== n) .&
                          (x.inv_on_hand.material .== p), :discarded] .+= onhand - max_inv
        end
    end
end

"""
    update_positions!(x::SupplyChainEnv, n::Int, p::Any)

Update inventory position and inventory level for material `p` at node `n`.
"""
function update_position!(x::SupplyChainEnv, n::Int, p::Any)
    making = sum(filter(j -> j.arc[end] == n && j.material == p, x.production).amount) #commited production order
    upstream = sum(filter(j -> j.period == x.period && j.arc[end] == n && j.material == p, x.inv_pipeline).level) #in-transit inventory
    onhand = filter(j -> j.period == x.period && j.node == n && j.material == p, x.inv_on_hand).level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.backlog
        if n in x.markets #find demand backlog
            backlog += filter(i -> i.period == x.period && i.node == n && i.material == p, x.demand).unfulfilled[1]
        else #find any unfulfilled replenishment request that was not reallocated
            backlog += sum(filter(i -> i.period == x.period && i.arc[1] == n && i.material == p && ismissing(i.reallocated), x.replenishments).unfulfilled)
        end
        backorder = sum(filter(i -> i.period == x.period && i.arc[end] == n && i.material == p && ismissing(i.reallocated), x.replenishments).unfulfilled)
    end
    push!(x.inv_level, [x.period, n, p, onhand - backlog]) #update inventory
    push!(x.inv_position, [x.period, n, p, onhand + making + upstream + backorder - backlog]) #update inventory
end

"""
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    for n in x.markets, p in x.materials
        dmnd_seq = get_prop(x.network, n, :demand_sequence)[p]
        freq = Bernoulli(get_prop(x.network, n, :demand_frequency)[p]) #demand frequency
        dmnd = get_prop(x.network, n, :demand_distribution)[p] #demand distribution
        if dmnd isa Sampleable
            dmnd = truncated(dmnd,0,Inf)
        end
        q = iszero(dmnd_seq) ? rand(freq) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
        demand = [x.period, n, p, q, 0., 0.] #initialize demand vector to store in df
        inv = filter(i -> i.period == x.period && i.node == n && i.material == p, x.inv_on_hand).level[1] #current inventory at node
        if x.backlog && x.period > 1 #add previous backlog to the quantity requested at the market
            q += filter(i -> i.period == x.period-1 && i.node == n && i.material == p, x.demand).unfulfilled[1]
        end
        demand[5] = min(q, inv) #sales
        demand[6] = max(q - inv, 0.) #unfilfilled
        push!(x.demand, demand) #update df

        #update end of period inventory (subtract sales)
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== n) .&
                      (x.inv_on_hand.material .== p), :level] .-= demand[5]

        #update inventory position
        update_position!(x, n, p)
    end
end

"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    for n in vertices(x.network)
        profit = 0 #initialize node profit
        for p in x.materials
            #get costs
            #holding cost
            hold_cost = get_prop(x.network, n, :holding_cost)[p]
            onhand = filter(j -> j.period == x.period && j.node == n && j.material == p, x.inv_on_hand).level[1]
            onhand = onhand == Inf ? 0 : onhand #avoid NaNs
            profit -= hold_cost * onhand
            #production cost
            if n in x.producers
                prod_cost = get_prop(x.network, n, :production_cost)[p]
                produced = sum(filter(j -> j.period == x.period && j.arc[1] == n && j.material == p, x.replenishments).accepted)
                profit -= prod_cost * produced
            #sales profit at markets (and penalize for unfulfilled demand)
            elseif n in x.markets
                sales_price = get_prop(x.network, n, :sales_price)[p]
                dmnd_penalty = get_prop(x.network, n, :demand_penalty)[p]
                sold, unfilled = filter(j -> j.period == x.period && j.node == n && j.material == p, x.demand)[1, [:sold, :unfulfilled]]
                profit += sales_price * sold
                profit -= dmnd_penalty * unfilled
            end
            #pay suppliers for received inventory and pay transportation cost
            for pred in inneighbors(x.network, n)
                price = get_prop(x.network, pred, n, :sales_price)[p]
                trans_cost = get_prop(x.network, pred, n, :transportation_cost)[p]
                #pay purchase of inventory
                purchased = filter(j -> j.arc == (pred, n) && j.material == p, arrivals).amount
                if !isempty(purchased)
                    profit -= purchased[1] * price
                end
                #pay transportation (pipeline hoding) cost (assume it is paid to a third party)
                intransit = filter(j -> j.period == x.period && j.arc == (pred, n) && j.material == p, x.inv_pipeline).level[1]
                profit -= intransit * trans_cost
            end
            #receive payment for delivered inventory
            for succ in outneighbors(x.network, n)
                price = get_prop(x.network, n, succ, :sales_price)[p]
                #receive payment for delivered inventory
                sold = filter(j -> j.arc == (n, succ) && j.material == p, arrivals).amount
                if !isempty(sold)
                    profit += sold[1] * price
                end
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end
