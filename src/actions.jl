"""
    (x::SupplyChainEnv)(action::Vector{T} where T <: Real)

Apply an `action` (replenishment requests) on the `SupplyChainEnv` and step
forward one simulation period.
"""
function (x::SupplyChainEnv)(action::Vector{T} where T <: Real)
    #validate action input
    # t0 = @elapsed begin
    @assert all(action .>= 0) "Reorder actions cannot be negative."
    @assert length(action) == length(x.materials)*ne(x.network) "Reorder action vector must have length num_products * num_edges."
    # end
    #increase period counter
    x.period += 1

    #intialize next period on-hand and pipeline inventories with previous inventories
    # t1 = @elapsed begin
    prev_inv_on_hand = filter([:period] => j -> j == x.period-1, x.inv_on_hand, view=true) #previous inventory level
    prev_ioh_group = groupby(prev_inv_on_hand, [:node, :material])
    prev_inv_pipeline = filter([:period] => j -> j == x.period-1, x.inv_pipeline, view=true) #previous inventory level
    prev_ip_group = groupby(prev_inv_pipeline, [:arc, :material])
    for p in x.materials
        for n in vertices(x.network)
            # prev = x.inv_on_hand |> @filter(i -> i.period == x.period-1 && i.node == n && i.material == p) |> DataFrame #previous inventory level
            # prev = filter([:period, :node, :material] => (j1, j2, j3) -> j1 == x.period-1 && j2 == n && j3 == p, x.inv_on_hand, view=true) #previous inventory level
            # push!(x.inv_on_hand, [x.period, n, p, prev[1,:level], 0]) #intialize with previous inventory levels
            push!(x.inv_on_hand, [x.period, n, p, prev_ioh_group[(node = n, material = p)][1,:level], 0]) #intialize with previous inventory levels
        end
        for a in edges(x.network)
            # prev = x.inv_pipeline |> @filter(i -> i.period == x.period-1 && i.arc == (a.src, a.dst) && i.material == p) |> DataFrame #previous inventory level
            # prev = filter([:period, :arc, :material] => (j1, j2, j3) -> j1 == x.period-1 && j2 == (a.src,a.dst) && j3 == p, x.inv_pipeline, view=true) #previous inventory level
            # push!(x.inv_pipeline, [x.period, (a.src,a.dst), p, prev[1,:level]]) #intialize with previous inventory levels
            push!(x.inv_pipeline, [x.period, (a.src,a.dst), p, prev_ip_group[(arc = (a.src,a.dst), material = p)][1,:level]]) #intialize with previous inventory levels
        end
    end
    # end
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
    # t5 = @elapsed begin
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)
    # end

    #update inventory positions for non-market nodes
    # t6 = @elapsed begin
    # for n in union(x.producers, x.distributors), p in x.materials
    #     update_position!(x, n, p)
    # end
    # end
    #markets open and demand occurs
    simulate_markets!(x)

    update_positions!(x)

    # println("$(x.period): $t0, $t1, $t2, $t3, $t4, $t5, $t6, $t7")

    if x.options[:evaluate_profit]
        #calculate profit at each node
        calculate_profit!(x, arrivals)

        #update reward (current profit). Should be revised for RL
        x.reward = sum(filter([:period] => j -> j == x.period, x.profit, view=true).value)
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
            push!(x.replenishments, [x.period, a, p, 0, 0, 0, 0, 0])
        end
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

    nonsources = [n for n in vertices(x.network) if !isempty(inneighbors(x.network, n))]
    supply_df = filter([:period] => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory at supplier
    supply_grp = groupby(supply_df, [:node, :material])
    for (i,p) in enumerate(mats) #loop by materials
        for n in nonsources #loop by nodes placing requests
            sup_priority = get_prop(x.network, n, :supplier_priority)[p] #get supplier priority list
            for (l, sup) in enumerate(sup_priority) #loop by supplier priority
                a = (sup, n) #arc
                j = findfirst(k -> k == a, arcs) #find index for that arc in the action matrix
                amount = act[i,j] #amount requested
                if x.options[:backlog] && x.period > 1 #add previous period's backlog
                    if x.options[:reallocate] && l == 1
                        amount += filter([:material, :arc, :period] => (k1, k2, k3) -> k1 == p && k2 == (sup_priority[end],n) && k3 == x.period-1, x.replenishments, view=true).unfulfilled[1]
                    elseif !x.options[:reallocate]
                        amount += filter([:material, :arc, :period] => (k1, k2, k3) -> k1 == p && k2 == a && k3 == x.period-1, x.replenishments, view=true).unfulfilled[1]
                    end
                end
                if iszero(amount) #continue to next iteration if request is 0
                    push!(x.replenishments, [x.period, a, p, requests[i,j], 0, 0, 0, 0])
                    continue 
                end
                # supply = filter([:period, :node, :material] => (k1, k2, k3) -> k1 == x.period && k2 == a[1] && k3 == p, x.inv_on_hand, view=true).level[1] #on_hand inventory at supplier
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
                        sup_pp = supply_grp[(node = a[1], material = mats[ii])].level[1]
                        # sup_pp = filter([:period, :node, :material] => (k1, k2, k3) -> k1 == x.period && k2 == a[1] && k3 == mats[ii], x.inv_on_hand, view=true).level[1] #supply of pp
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

                #chekc if some was not accepted and reallocate
                # unfulfilled = amount - accepted #amount not accepted
                new_alloc = missing
                if amount > 0 #reallocate unfulfilled request to next priority supplier
                    # @warn "Replenishment request made in period $(x.period) by node $(a[2]) to node $(a[1]) for material $p was reduced from $amount to $accepted due to insufficient production capacity or insufficient inventory."
                    if x.options[:reallocate]
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
    produced = filter([:lead] => i -> i <= 0, x.production, view=true) #find active production that has completed
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
    filter!([:lead] => i -> i > 0, x.production) #remove produced inventory that was shipped (and any zero values)
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    arrivals = filter([:lead] => i -> i <= 0, x.shipments, view=true) #find active shipments with 0 lead time
    for i in 1:nrow(arrivals)
        a, p, amount = arrivals[i, [:arc, :material, :amount]]
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== a[end]) .&
                      (x.inv_on_hand.material .== p), :level] .+= amount
        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                       (x.inv_pipeline.material .== p), :level] .-= amount
    end
    filter!([:lead] => i -> i > 0, x.shipments) #remove shipments that arrived (and any zero values)
    return arrivals
end

"""
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    supply_df = filter([:period] => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at supplier
    supply_grp = groupby(supply_df, [:node, :material])
    for n in vertices(x.network)
        node_max_inv = get_prop(x.network, n, :inventory_capacity)
        for p in x.materials
            max_inv = node_max_inv[p]
            isinf(max_inv) && continue
            onhand = supply_grp[(node = n, material = p)]
            # onhand = filter([:period, :node, :material] => (j1, j2, j3) -> j1 == x.period && j2 == n && j3 == p, x.inv_on_hand, view=true).level[1] #on_hand inventory
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
end

"""
    update_positions!(x::SupplyChainEnv, n::Int, p::Any)

Update inventory position and inventory level for material `p` at node `n`.
"""
function update_position!(x::SupplyChainEnv, n::Int, p::Any)
    making = sum(filter([:arc, :material] => (j1, j2) -> j1[end] == n && j2 == p, x.production, view=true).amount) #commited production order
    upstream = sum(filter([:period, :arc, :material] => (j1, j2, j3) -> j1 == x.period && j2[end] == n && j3 == p, x.inv_pipeline, view=true).level) #in-transit inventory
    onhand = filter([:period, :node, :material] => (j1, j2, j3) -> j1 == x.period && j2 == n && j3 == p, x.inv_on_hand, view=true).level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        if n in x.markets #find demand backlog
            backlog += filter([:period, :node, :material] => (i1, i2, i3) -> i1 == x.period && i2 == n && i3 == p, x.demand, view=true).unfulfilled[1]
        else #find any unfulfilled replenishment request that was not reallocated
            backlog += sum(filter([:period, :arc, :material, :reallocated] => (i1, i2, i3, i4) -> i1 == x.period && i2[1] == n && i3 == p && ismissing(i4), x.replenishments, view=true).unfulfilled)
        end
        backorder = sum(filter([:period, :arc, :material, :reallocated] => (i1, i2, i3, i4) -> i1 == x.period && i2[end] == n && i3 == p && ismissing(i4), x.replenishments, view=true).unfulfilled)
    end
    push!(x.inv_level, [x.period, n, p, onhand - backlog]) #update inventory
    push!(x.inv_position, [x.period, n, p, onhand + making + upstream + backorder - backlog]) #update inventory
end


"""
    update_positions!(x::SupplyChainEnv)

Update inventory position and inventory level for all materials all nodes.
"""
function update_positions!(x::SupplyChainEnv)
    supply_df = filter([:period] => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at supplier
    supply_grp = groupby(supply_df, [:node, :material])
    pipeline_df = filter([:period] => j -> j == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, :material)
    dmnd_df = filter([:period] => i -> i == x.period, x.demand, view=true) #demand at markets
    dmnd_grp = groupby(dmnd_df, [:node, :material])
    orders_df = filter([:period, :reallocated] => (i1, i2) -> i1 == x.period && ismissing(i2), x.replenishments, view=true)
    orders_grp = groupby(orders_df, :material)

    for n in vertices(x.network), p in x.materials
        making = sum(filter([:arc, :material] => (j1, j2) -> j1[end] == n && j2 == p, x.production, view=true).amount) #commited production order
        upstream = sum(filter([:arc] => j -> j[end] == n, pipeline_grp[(material = p,)], view=true).level) #in-transit inventory
        # upstream = sum(filter([:period, :arc, :material] => (j1, j2, j3) -> j1 == x.period && j2[end] == n && j3 == p, x.inv_pipeline, view=true).level) #in-transit inventory
        onhand = supply_grp[(node = n, material = p)].level[1] #on_hand inventory
        # onhand = filter([:period, :node, :material] => (j1, j2, j3) -> j1 == x.period && j2 == n && j3 == p, x.inv_on_hand, view=true).level[1] #on_hand inventory
        backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
        backlog = 0 #initialize backlog for orders placed by successors
        if x.options[:backlog]
            if n in x.markets #find demand backlog
                backlog += dmnd_grp[(node = n, material = p)].unfulfilled[1]
                # backlog += filter([:period, :node, :material] => (i1, i2, i3) -> i1 == x.period && i2 == n && i3 == p, x.demand, view=true).unfulfilled[1]
            else #find any unfulfilled replenishment request that was not reallocated
                backlog += sum(filter([:arc] => i -> i[1] == n, orders_grp[(material = p,)], view=true).unfulfilled)
                # backlog += sum(filter([:period, :arc, :material, :reallocated] => (i1, i2, i3, i4) -> i1 == x.period && i2[1] == n && i3 == p && ismissing(i4), x.replenishments, view=true).unfulfilled)
            end
            backorder = sum(filter([:arc] => i -> i[end] == n, orders_grp[(material = p,)], view=true).unfulfilled)
            # backorder = sum(filter([:period, :arc, :material, :reallocated] => (i1, i2, i3, i4) -> i1 == x.period && i2[end] == n && i3 == p && ismissing(i4), x.replenishments, view=true).unfulfilled)
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
    supply_df = filter([:period] => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at supplier
    supply_grp = groupby(supply_df, [:node, :material])
    last_dmnd = filter([:period] => i -> i == x.period-1, x.demand, view=true)
    last_d_group = groupby(last_dmnd, [:node, :material])
    for n in x.markets, p in x.materials
        dmnd_seq = get_prop(x.network, n, :demand_sequence)[p]
        freq = Bernoulli(get_prop(x.network, n, :demand_frequency)[p]) #demand frequency
        dmnd = get_prop(x.network, n, :demand_distribution)[p] #demand distribution
        q = iszero(dmnd_seq) ? rand(freq) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
        demand = [x.period, n, p, q, 0., 0.] #initialize demand vector to store in df
        if x.options[:backlog] && x.period > 1 #add previous backlog to the quantity requested at the market
            q += last_d_group[(node = n, material = p)].unfulfilled[1]
            # q += filter([:period, :node, :material] => (i1, i2, i3) -> i1 == x.period-1 && i2 == n && i3 == p, x.demand, view=true).unfulfilled[1]
        end
        if q > 0
            inv = supply_grp[(node = n, material = p)].level[1]
            # inv = filter([:period, :node, :material] => (i1, i2, i3) -> i1 == x.period && i2 == n && i3 == p, x.inv_on_hand, view=true).level[1] #current inventory at node
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

        # #update inventory position
        # update_position!(x, n, p)
    end
end

"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    supply_df = filter([:period, :level] => (j1, j2) -> j == x.period && !isinf(j2), x.inv_on_hand, view=true) #on_hand inventory at supplier
    supply_grp = groupby(supply_df, [:node, :material])
    orders_df = filter([:period] => i -> i == x.period, x.replenishments, view=true)
    orders_grp = groupby(orders_df, :material)
    sales_df = filter([:period] => j -> j == x.period, x.demand, view=true)
    sales_grp = groupby(sales_df, [:node, :material])
    pipeline_df = filter([:period] => j -> j == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])
    for n in vertices(x.network)
        profit = 0 #initialize node profit
        for p in x.materials
            #get costs
            #holding cost
            hold_cost = get_prop(x.network, n, :holding_cost)[p]
            if hold_cost > 0
                onhand = supply_grp[(node = n, material = p)].level
                # onhand = filter([:period, :node, :material, :level] => (j1, j2, j3, j4) -> j1 == x.period && j2 == n && j3 == p && !isinf(j4), x.inv_on_hand, view=true).level
                if !isempty(onhand)
                    profit -= hold_cost * onhand[1]
                end
            end
            #production cost
            if n in x.producers
                prod_cost = get_prop(x.network, n, :production_cost)[p]
                if prod_cost > 0
                    produced = sum(filter([:arc] => j -> j[1] == n, orders_grp[(material = p,)], view=true).accepted)
                    # produced = sum(filter([:period, :arc, :material] => (j1, j2, j3) -> j1 == x.period && j2[1] == n && j3 == p, x.replenishments, view=true).accepted)
                    profit -= prod_cost * produced
                end
            #sales profit at markets (and penalize for unfulfilled demand)
            elseif n in x.markets
                sales_price = get_prop(x.network, n, :sales_price)[p]
                dmnd_penalty = get_prop(x.network, n, :demand_penalty)[p]
                if sales_price > 0 || dmnd_penalty > 0
                    sold, unfilled = sales_grp[(node = n, material = p)][1, [:sold, :unfulfilled]]
                    # sold, unfilled = filter([:period, :node, :material] => (j1, j2, j3) -> j1 == x.period && j2 == n && j3 == p, x.demand, view=true)[1, [:sold, :unfulfilled]]
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
                    # intransit = filter([:period, :arc, :material] => (j1, j2, j3) -> j1 == x.period && j2 == (pred, n) && j3 == p, x.inv_pipeline, view=true).level[1]
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
