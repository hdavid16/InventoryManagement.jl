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
    initialize_inventories!(x)
    
    #move active shipments forward one period
    x.shipments.lead .-= 1

    # #move active production forward one period
    # x.production.lead .-= 1

    #reshape action
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    act = reshape(action, (length(x.materials), length(arcs)))

    #place requests
    place_requests!(x, act, arcs)

    # #update production and associated pipeline inventories
    # update_production!(x)

    #update on hand and pipeline inventories due to arrived shipments
    arrivals = update_shipments!(x)

    #discard any excess inventory
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)

    #markets open and demand occurs
    simulate_markets!(x)

    #update inventories at each node and echelon
    update_inventories!(x)

    if x.options[:evaluate_profit]
        #calculate profit at each node
        calculate_profit!(x, arrivals)
        #update reward (current profit). Should be revised for RL
        x.reward = sum(filter(:period => j -> j == x.period, x.profit, view=true).value)
    end
end

"""
    initialize_inventories!(x::SupplyChainEnv)

Initialize on-hand and pipeline inventories with those from the previous period.
"""
function initialize_inventories!(x::SupplyChainEnv)
    #find previous period's inventories
    prev_inv_on_hand = filter(:period => j -> j == x.period-1, x.inv_on_hand, view=true) #previous inventory level
    prev_onhand_grp = groupby(prev_inv_on_hand, [:node, :material])
    prev_inv_pipeline = filter(:period => j -> j == x.period-1, x.inv_pipeline, view=true) #previous inventory level
    prev_pipeln_grp = groupby(prev_inv_pipeline, [:arc, :material])

    for mat in x.materials
        for n in vertices(x.network)
            prev = prev_onhand_grp[(node = n, material = mat)][1,:level] #previous on hand inventory
            push!(x.inv_on_hand, [x.period, n, mat, prev, 0]) #intialize with previous inventory levels
        end
        for a in edges(x.network)
            prev = prev_pipeln_grp[(arc = (a.src,a.dst), material = mat)][1,:level] #previous pipeline inventory
            push!(x.inv_pipeline, [x.period, (a.src,a.dst), mat, prev]) #intialize with previous inventory levels
        end
    end
end

"""
    place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)

Place inventory replenishment requests throughout the network.
"""
function place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)
    #get previous orders
    last_orders_grp = previous_orders(x)
    #exit if no action
    if iszero(act)
        exit_order(x, arcs, last_orders_grp)
        return
    end

    #extract info
    mats = x.materials
    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    #sample lead times
    leads = Dict((a,mat) => rand(get_prop(x.network, a, :lead_time)[mat]) for a in edges(x.network), mat in mats)
    #non source nodes
    nonsources = [n for n in vertices(x.network) if !isempty(inneighbors(x.network, n))]
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #place requests
    for (i,mat) in enumerate(mats) #loop by materials
        for req in nonsources #loop by nodes placing requests
            sup_priority = get_prop(x.network, req, :supplier_priority)[mat] #get supplier priority list
            for sup in sup_priority #loop by supplier priority
                a = (sup, req) #arc
                j = findfirst(k -> k == a, arcs) #find index for that arc in the action matrix
                lead = float(leads[Edge(a...), mat]) #sampled lead time
                
                #amount requested
                amount = act[i,j]
                #add previous period's backlog
                amount = add_backlog!(x, a..., mat, amount, sup_priority, last_orders_grp)

                #continue to next iteration if request is 0
                if iszero(amount) 
                    push!(x.replenishments, [x.period, a, mat, 0, 0, 0, 0, missing])
                    continue 
                end

                #accept or adjust requests
                if sup in x.producers #supplier is a plant
                    amount, accepted = producer_fulfillment!(x, a..., mat, amount, lead, supply_grp, pipeline_grp, capacities)
                else #supplier is a distributor
                    amount, accepted = fulfill_from_stock!(x, a..., mat, amount, lead, supply_grp, pipeline_grp)
                end

                #check if some was not accepted and reallocate
                new_alloc = reallocate_request!(x, a..., amount, sup_priority, arcs, act)

                #store shipments and lead times
                if accepted > 0 #update active shipments
                    if amount > 0
                        x.replenishments[end, :unfulfilled] = amount
                        x.replenishments[end, :reallocated] = new_alloc
                    end
                else #no request made
                    push!(x.replenishments, [x.period, a, mat, amount, 0, 0, amount, new_alloc])
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
    previous_orders(x::SupplyChainEnv)

Returns a view of a grouped dataframe of the replenishment orders for the previous period, filtered by arc and material.
"""
function previous_orders(x::SupplyChainEnv)
    last_orders_grp = missing
    if x.options[:backlog] && x.period > 1
        last_orders_df = filter(:period => i -> i == x.period-1, x.replenishments, view=true) #replenishment orders from previous period
        last_orders_grp = groupby(last_orders_df, [:arc, :material]) #group by material and arc
    end

    return last_orders_grp
end

"""
    exit_order(x, arcs, last_orders_grp)

Abort order placement.
"""
function exit_order(x::SupplyChainEnv, arcs::Vector, last_orders_grp::Union{Missing, GroupedDataFrame})
    backlog = 0
    for a in arcs, mat in x.materials
        if x.options[:backlog] && x.period > 1
            sup_priority = get_prop(x.network, a[2], :supplier_priority)[mat] #get supplier priority list
            if x.options[:reallocate] && a[1] == sup_priority[1]
                backlog = last_orders_grp[(arc = (sup_priority[end],a[2]), material = mat)].unfulfilled[end]
            elseif !x.options[:reallocate]
                backlog = last_orders_grp[(arc = a, material = mat)].unfulfilled[end]
            end
        end
        push!(x.replenishments, [x.period, a, mat, 0, 0, 0, backlog, missing])
    end
end

"""
    add_backlog!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, amount::Number, sup_priority::Vector, last_orders_grp::Union{Missing, GroupedDataFrame})

Add backlog to requested amount.
"""
function add_backlog!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, amount::Float64, sup_priority::Vector, last_orders_grp::Union{Missing, GroupedDataFrame})
    l = findfirst(i -> i == sup, sup_priority) #location of sup on priority ranking
    #add previous period's backlog
    if x.options[:backlog] && x.period > 1 
        if x.options[:reallocate] && l == 1 #unfulfilled from lowest priority supplier does not get reallocated explicitly (so add it here)
            amount += last_orders_grp[(arc = (sup_priority[end], req), material = mat)].unfulfilled[1]
        elseif !x.options[:reallocate] #add the previous backlog at the node since reallocaiton was not done
            amount += last_orders_grp[(arc = (sup, req), material = mat)].unfulfilled[1]
        end
    end

    return amount
end

"""
    producer_fulfillment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict)

Order fulfillment by plant node.
"""
function producer_fulfillment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict)
    #first, try to fulfill with on-hand inventory
    amount, accepted_inv = fulfill_from_stock!(x, src, dst, mat, amount, lead, supply_grp, pipeline_grp)
    # #next, try to satisfy with by-products
    # amount, accepted_sched = fulfill_from_coproduction!(x, src, dst, mat, amount, lead)
    #finally, try to satisfy with material production
    amount, accepted_prod = fulfill_from_production!(x, src, dst, mat, amount, lead, capacities, supply_grp, pipeline_grp)

    #total accepted request
    # accepted = accepted_inv + accepted_sched + accepted_prod
    accepted = accepted_inv + accepted_prod
    return amount, accepted
end

"""
    fulfill_from_stock!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)

Fulfill request from on-hand inventory.
"""
function fulfill_from_stock!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)
    #available supply
    supply = supply_grp[(node = src, material = mat)].level[1]
    
    #try to satisfy with on hand inventory first
    original_amount = amount #NOTE: change from v0.3.4 = The original amount will now include any backlog (affects service measures), and get's updated for each fulfillment type (stock, production)
    accepted_inv = min(amount, supply) 
    if accepted_inv > 0
        amount -= accepted_inv #get new amount pending
        supply_grp[(node = src, material = mat)].level[1] -= accepted_inv #remove inventory from site
        #ship material
        ship!(x, src, dst, mat, accepted_inv, original_amount, lead, pipeline_grp)
    end

    return amount, accepted_inv
end

# """
#     fulfill_from_coproduction!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64)

# Fulfill request from scheduled coproduction.
# """
# function fulfill_from_coproduction!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64)
#     sort!(x.production, :lead) #sort to use scheduled supply that is closest to finishing
#     sched_supply = filter([:arc, :material] => (k1, k2) -> k1 == (src,src) && k2 == mat, x.production, view=true) #check if any inventory is scheduled for production, but not commited to downstream node
#     accepted_sched_k = [] #try to satisfy with scheduled production that is not commited
#     if !isempty(sched_supply)
#         for (kx, k) in enumerate(sched_supply.amount)
#             original_amount = amount #save amount requested
#             accepted_k = min(amount, k) #satisfy if possible with scheduled supply
#             push!(accepted_sched_k, accepted_k) #store accepted amount
#             push!(x.production, ((src,dst), mat, accepted_k, sched_supply[kx, :lead])) #add reallocatied production to x.production
#             if accepted_k > 0
#                 amount -= accepted_k #update pending `amount`
#                 sched_supply[kx, :amount] -= accepted_k #update x.production amount (reduce by reallocated amount)
#                 sched_lead = sched_supply[kx, :lead]
#                 ship!(x, src, dst, mat, accepted_k, original_amount, lead, sched_lead)
#             end
#         end
#     end
#     accepted_sched = reduce(+, accepted_sched_k, init = 0)

#     return amount, accepted_sched
# end

"""
    fulfill_from_production!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, capacities::Dict, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, amount::Float64, lead::Float64, capacities::Dict, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)
    #extract info
    mats = x.materials
    bom = get_prop(x.network, src, :bill_of_materials)
    # bom = x.bill_of_materials
    # i = findfirst(j -> j == mat, mats)

    #commit production at plant
    capacity = [] #get production capacity
    mat_supply = [] #store max capacity based on raw material consumption for each raw material
    rmats = findall(k -> k < 0, bom[:,mat]) #indices of raw materials involved with production of mat
    cmats = findall(k -> k > 0, bom[:,mat]) #indices for co-products
    rmat_names = names(bom,1)[rmats] #names of raw materials
    cmat_names = names(bom,1)[cmats] #names of co-products
    if !isempty(rmats) #only add capacity if there is a non-zero bom for that material
        push!(capacity, capacities[src][mat])
    else
        push!(capacity, 0)
    end
    for rmat in rmat_names
        sup_pp = supply_grp[(node = src, material = rmat)].level[1] #supply of material involved in BOM
        push!(mat_supply, - sup_pp / bom[rmat,mat]) #only account for raw materials that are in the BOM
    end 
    for cmat in cmat_names #add capacity constraint for any co-products (scaled by stoichiometry)
        push!(capacity, capacities[src][rmat] / bom[cmat,mat])
    end
    accepted_prod = min(amount, capacity..., mat_supply...) #fulfill remaining with available capacity

    #schedule production & shipment and update inventories
    if accepted_prod > 0
        # prod_time = get_prop(x.network, src, :production_time)[mat] #production time
        original_amount = amount
        amount -= accepted_prod #update pending amount
        ship!(x, src, dst, mat, accepted_prod, original_amount, lead, pipeline_grp) #schedule shipment
        # push!(x.production, ((src,dst), mat, accepted_prod, prod_time))
        capacities[src][mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
        for rmat in rmat_names #reactant consumed
            supply_grp[(node = src, material = rmat)].level[1] += accepted_prod * bom[rmat,mat]
        end
        for cmat in cmat_names #coproducts scheduled for production
            # push!(x.production, ((src,src), mats[ii], accepted_prod*bom[ii,i], prod_time))
            coproduction = accepted_prod*bom[cmat,mat]
            ship!(x, src, dst, cmat, coproduction, missing, lead, pipeline_grp) #schedule shipment of co-product (original production is missing since this is a coproduct)
            capacities[src][cmat] -= coproduction #update coproduct production capacity
        end
        # #schedule shipment
        # ship!(x, src, dst, mat, accepted_prod, original_amount, lead)#, float(prod_time))
    end

    return amount, accepted_prod
end

"""
    ship!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, original_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame, delay::Float64 = 0.)

Schedule shipment of material.
"""
function ship!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, original_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame, delay::Float64 = 0.)
    push!(x.shipments, [(src,dst), mat, accepted_amount, lead + delay])
    push!(x.replenishments, [x.period, (src,dst), mat, original_amount, accepted_amount, lead + delay, 0, missing])
    pipeline_grp[(arc = (src, dst), material = mat)].level[1] += accepted_amount #add inventory to pipeline
end

"""
    reallocate_request!(x::SupplyChainEnv, sup::Int, req::Int, amount::Float64, sup_priority::Vector, arcs::Vector, act::Array)

Reallocate order to next supplier in the priority list.
"""
function reallocate_request!(x::SupplyChainEnv, sup::Int, req::Int, amount::Float64, sup_priority::Vector, arcs::Vector, act::Array)
    new_alloc = missing
    if amount > 0 && x.options[:reallocate] #reallocate unfulfilled request to next priority supplier
        next_sup = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
        if next_sup <= length(sup_priority) #check that there is a next one in line
            new_sup = sup_priority[next_sup] #next supplier in line
            new_alloc = (new_sup, req) #store new arc where reallocated
            jj = findfirst(k -> k == new_alloc, arcs) #find index for that arc in the action matrix
            act[i,jj] += amount #add unfulfilled to the next supplier in the line
        end
    end

    return new_alloc
end

# """
#     update_production!(x::SupplyChainEnv)

# Update completed production at each producer node and send to the arc that it
#     was commited to.
# """
# function update_production!(x::SupplyChainEnv)
#     #extract info
#     mats = x.materials
#     bom = x.bill_of_materials

#     #filter data
#     #get on hand inventory
#     supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
#     supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
#     #get pipeline inventory
#     pipeline_df = filter(:period => k -> k == x.period, x.inv_pipeline, view=true)
#     pipeline_grp = groupby(pipeline_df, [:arc, :material])

#     #find production that has completed
#     produced = filter([:lead, :amount] => (i1, i2) -> i1 <= 0 && i2 > 0, x.production) #find active production that has completed
#     for i in 1:nrow(produced) #add produced material to pipeline
#         a, mat, amount = produced[i, [:arc, :material, :amount]]
#         #restore production capacities
#         capacity = get_prop(x.network, a[1], :production_capacity)
#         capacity[mat] += amount
#         j = findfirst(k -> k==mat, mats)
#         imats = findall(k -> k > 0, bom[:,j])
#         for ii in imats
#             capacity[mats[ii]] += amount*bom[ii,j]
#         end
#         set_prop!(x.network, a[1], :production_capacity, capacity)
#         #update inventories
#         if a[1] != a[2] #send produced material down pipeline
#             pipeline_grp[(arc = a, material = mat)].level[1] += amount
#         else #send produced material to storage at that node
#             supply_grp[(node = a[1], material = mat)].level[1] += amount
#         end
#     end
#     filter!([:lead, :amount] => (i1,i2) -> i1 > 0 && i2 > 0, x.production) #remove produced inventory that was shipped (and any zero values)
# end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    #filter data
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #find active shipments with 0 lead time
    arrivals = filter(:lead => i -> i <= 0, x.shipments) 
    for i in 1:nrow(arrivals)
        a, mat, amount = arrivals[i, [:arc, :material, :amount]]
        supply_grp[(node = a[end], material = mat)].level[1] += amount
        pipeline_grp[(arc = a, material = mat)].level[1] -= amount
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
        for mat in x.materials
            max_inv = node_max_inv[mat]
            (iszero(max_inv) || isinf(max_inv)) && continue #skip iteration if node doesn't store that material or has infinite capacity
            onhand = onhand_grp[(node = n, material = mat)].level[1] #onhand inventory
            if onhand > max_inv
                onhand_grp[(node = n, material = mat)].level[1] = max_inv
                onhand_grp[(node = n, material = mat)].discarded[1] += onhand - max_inv
            end
        end
    end
end


"""
    update_inventories!(x::SupplyChainEnv)

Update inventory position and inventory level for all materials and nodes.
"""
function update_inventories!(x::SupplyChainEnv)
    #filter data
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, :material) #group by material
    dmnd_df = filter(:period => i -> i == x.period, x.demand, view=true) #demand at markets
    dmnd_grp = groupby(dmnd_df, [:node, :material]) #group by node and material
    orders_df = filter([:period, :reallocated] => (i1, i2) -> i1 == x.period && ismissing(i2), x.replenishments, view=true) #replenishment orders from current period
    orders_grp = groupby(orders_df, :material) #group by material

    #get material conversion
    # dmnd_quant = combine(groupby(dmnd_df, :material), :sold => sum) #add up demand for each material
    # dmnd_dict = Dict(dmnd_quant.material .=> dmnd_quant.sold_sum) #convert to dictionary
    # conversion_dict = get_prop(x.network, :conversion_dictionary)

    #initialize echelon inventory positions
    initialize_echelons!(x)
    ech_df = filter(:period => j -> j == x.period, x.ech_position, view=true) #echelon position at each node
    ech_grp = groupby(ech_df, [:node, :material]) #group by material
    
    #loop through nodes and update inventory levels, positions, and echelons
    for n in vertices(x.network), mat in x.materials
        ilevel, ipos0, ipos = inventory_components(x, n, mat, pipeline_grp, onhand_grp, dmnd_grp, orders_grp)
        #update inventory
        push!(x.inv_level, [x.period, n, mat, ilevel]) 
        #update inventory
        push!(x.inv_position, [x.period, n, mat, ipos]) 
        #update echelons
        update_echelons!(x, n, mat, ipos, ipos0, ech_grp)#, dmnd_dict, conversion_dict)
    end

end

"""
    initialize_echelons!(x::SupplyChainEnv)

Initialize echelon positions at the current period to 0.
"""
function initialize_echelons!(x::SupplyChainEnv)
    for n in vertices(x.network), mat in x.materials
        if get_prop(x.network, n, :inventory_capacity)[mat] > 0
            push!(x.ech_position, [x.period, n, mat, 0])
        end
    end
end

"""
    inventory_components(x::SupplyChainEnv, n::Int, mat::Symbol, pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, dmnd_grp::GroupedDataFrame, orders_grp::GroupedDataFrame)

Extract components to determine inventory paramters.
"""
function inventory_components(x::SupplyChainEnv, n::Int, mat::Symbol, pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, dmnd_grp::GroupedDataFrame, orders_grp::GroupedDataFrame)
    # making = reduce(+,filter([:arc, :material] => (j1, j2) -> j1[end] == n && j2 == mat, x.production, view=true).amount, init=0) #commited production order
    upstream = sum(filter(:arc => j -> j[end] == n, pipeline_grp[(material = mat,)], view=true).level) #in-transit inventory
    onhand = onhand_grp[(node = n, material = mat)].level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        if n in x.markets #find demand backlog
            backlog += dmnd_grp[(node = n, material = mat)].unfulfilled[1]
        else #find any unfulfilled replenishment request that was not reallocated
            backlog += sum(filter(:arc => i -> i[1] == n, orders_grp[(material = mat,)], view=true).unfulfilled)
        end
        backorder = sum(filter(:arc => i -> i[end] == n, orders_grp[(material = mat,)], view=true).unfulfilled)
    end
    ilevel = onhand - backlog #inventory level
    ipos0 = onhand + upstream + backorder #inventory position without backlog
    # ipos0 = onhand + making + upstream + backorder #inventory position without backlog
    ipos = ipos0 - backlog #include backlog in inventory position

    return ilevel, ipos0, ipos
end

"""
    update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ipos0::Float64, ech_grp::GroupedDataFrame)

Update echelon positions for current time period.
"""
function update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ipos0::Float64, ech_grp::GroupedDataFrame)#, dmnd_dict::Dict, conversion_dict::Dict)
    #identify which echelons have been affected and add to these
    for ech in findall(i -> n in i, x.echelons)
        if get_prop(x.network, ech, :inventory_capacity)[mat] > 0 #only add to echelon if that node holds that material
            if n in x.markets #backlog is only added for market nodes (to avoid double counting with backorder)
                ech_grp[(node = ech, material = mat)].level[1] += ipos
            else
                ech_grp[(node = ech, material = mat)].level[1] += ipos0
            end
        end
    end
    # #adjust all intermediates/raw materials for demand at the markets to correct the echelon position to account for these
    # for f in setdiff(x.products, [mat])
    #     if get_prop(x.network, n, :inventory_capacity)[mat] > 0
    #         ech_grp[(node = n, material = mat)].level[1] += dmnd_dict[f] * conversion_dict[mat,f] #subtract "equivalent" raw material sold
    #     end
    # end
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
    for n in x.markets, mat in x.materials
        dmnd_seq = get_prop(x.network, n, :demand_sequence)[mat]
        dprob = Bernoulli(1/get_prop(x.network, n, :demand_frequency)[mat]) #demand probability (probability of ordering is 1/demand period; aka, once every x days)
        dmnd = get_prop(x.network, n, :demand_distribution)[mat] #demand distribution
        q = iszero(dmnd_seq) ? rand(dprob) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
        if x.options[:backlog] && x.period > 1 #add previous backlog to the quantity requested at the market
            q += last_d_group[(node = n, material = mat)].unfulfilled[1]
        end
        inv = onhand_grp[(node = n, material = mat)].level[1]
        sold = min(q, inv) #sales
        unfilled = max(q - inv, 0.) #unfilfilled
        push!(x.demand, [x.period, n, mat, q, sold, unfilled]) #update df

        #update end of period inventory (subtract sales)
        onhand_grp[(node = n, material = mat)].level[1] -= sold
    end
end

"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    #filter data
    on_hand_df = filter([:period, :level] => (j1, j2) -> j1 == x.period && !isinf(j2), x.inv_on_hand, view=true) #on_hand inventory 
    onhand_grp = groupby(on_hand_df, [:node, :material])
    orders_df = filter(:period => i -> i == x.period, x.replenishments, view=true) #replenishment orders
    orders_grp = groupby(orders_df, :material)
    sales_df = filter(:period => j -> j == x.period, x.demand, view=true) #sales at markets
    sales_grp = groupby(sales_df, [:node, :material])
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #evaluate node profit
    for n in vertices(x.network)
        profit = 0. #initialize node profit
        for mat in x.materials
            #get costs
            #holding cost
            holding_cost = holding_costs(x, n, mat, onhand_grp)
            #production cost
            # if n in x.producers
            #     profit += production_costs(x, n, mat, orders_grp)
            #sales profit at markets (and penalize for unfulfilled demand)
            if n in x.markets
                profit += sum(sales_and_penalties(x, n, mat, sales_grp))
            end
            #pay suppliers for received inventory and pay transportation/production cost
            for pred in inneighbors(x.network, n)
                profit += accounts_payable(x, pred, n, mat, arrivals, pipeline_grp)
            end
            #receive payment for delivered inventory
            for succ in outneighbors(x.network, n)
                profit += accounts_receivable(x, n, succ, mat, arrivals)
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end

"""
    holding_costs(x::SupplyChainEnv, n::Int, mat::Symbol, onhand_grp::GroupedDataFrame)

Calculate holding costs (negative).
"""
function holding_costs(x::SupplyChainEnv, n::Int, mat::Symbol, onhand_grp::GroupedDataFrame)
    #holding cost
    hold_cost = get_prop(x.network, n, :holding_cost)[mat]
    c = 0
    if hold_cost > 0
        onhand = onhand_grp[(node = n, material = mat)].level
        if !isempty(onhand)
            c -= hold_cost * onhand[1]
        end
    end

    return c
end

# """
#     production_costs(x::SupplyChainEnv, n::Int, mat::Symbol, orders_grp::GroupedDataFrame)

# Calculate produciton costs (negative).
# """
# function production_costs(x::SupplyChainEnv, n::Int, mat::Symbol, orders_grp::GroupedDataFrame)
#     prod_cost = get_prop(x.network, n, :production_cost)[mat]
#     c = 0
#     if prod_cost > 0
#         produced = sum(filter([:arc] => j -> j[1] == n, orders_grp[(material = mat,)], view=true).accepted)
#         c -= prod_cost * produced
#     end

#     return c
# end

"""
    sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Symbol, sales_grp::GroupedDataFrame)

Calculate sales (positive) and unfulfilled demand penalties (negative)
"""
function sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Symbol, sales_grp::GroupedDataFrame)
    sales_price = get_prop(x.network, n, :sales_price)[mat]
    dmnd_penalty = get_prop(x.network, n, :demand_penalty)[mat]
    s = 0 #sales
    p = 0 #penalties
    if sales_price > 0 || dmnd_penalty > 0
        sold, unfilled = sales_grp[(node = n, material = mat)][1, [:sold, :unfulfilled]]
        s += sales_price * sold
        p -= dmnd_penalty * unfilled
    end

    return s, p
end

"""
    accounts_payable(x::SupplyChainEnv, sup::Int, n::Int, mat::Symbol, arrivals::DataFrame, pipeline_grp::GroupedDataFrame)

Calculate payments to suppliers and transportation costs (negative).
"""
function accounts_payable(x::SupplyChainEnv, sup::Int, n::Int, mat::Symbol, arrivals::DataFrame, pipeline_grp::GroupedDataFrame)
    price = get_prop(x.network, sup, n, :sales_price)[mat]
    trans_cost = get_prop(x.network, sup, n, :transportation_cost)[mat]
    pipe_holding_cost = get_prop(x.network, sup, n, :pipeline_holding_cost)[mat]
    ap = 0 #accounts payable
    if price > 0 || trans_cost > 0 #pay purchase of inventory and transportation cost (assume it is paid to a third party)
        purchased = filter([:arc, :material] => (j1, j2) -> j1 == (sup, n) && j2 == mat, arrivals, view=true).amount
        if !isempty(purchased)
            ap -= purchased[1] * price
            ap -= purchased[1] * trans_cost
        end
    end
    if pipe_holding_cost > 0 #pay pipeline holding cost (paid for in-transit inventory in the current period)
        intransit = pipeline_grp[(arc = (sup, n), material = mat)].level[1]
        ap -= intransit * pipe_holding_cost
    end

    return ap
end

"""
    accounts_receivable(x::SupplyChainEnv, n::Int, req::Int, mat::Symbol, arrivals::DataFrame)

Calculate payment for sales to downstream nodes (positive).
"""
function accounts_receivable(x::SupplyChainEnv, n::Int, req::Int, mat::Symbol, arrivals::DataFrame)
    price = get_prop(x.network, n, req, :sales_price)[mat]
    ar = 0
    if price > 0 #receive payment for delivered inventory
        sold = filter([:arc, :material] => (j1, j2) -> j1 == (n, req) && j2 == mat, arrivals, view=true).amount
        if !isempty(sold)
            ar += sold[1] * price
        end
    end

    return ar
end