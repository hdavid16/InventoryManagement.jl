"""
    (x::SupplyChainEnv, action::Vector{T} where T <: Real)

Convert a replenishment order vector into a NamedArray indicating how much of
each material (rows) is being requested on each arc (columns).
"""
function show_action(x::SupplyChainEnv, action::Vector{T} where T <: Real)
    mats = x.materials
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    action_named_array = NamedArray(
        reshape(action, (length(mats), length(arcs))),
        (x.materials, arcs),
        (:material, :arc)
    )
    
    return action_named_array
end

"""
    (x::SupplyChainEnv)(action::Vector{T} where T <: Real)

Apply an `action` (replenishment requests) on the `SupplyChainEnv` and step
forward one simulation period.
"""
function (x::SupplyChainEnv)(action::Vector{T} where T <: Real)
    #validate action input and reshape to matrix
    @assert all(action .>= 0) "Reorder actions cannot be negative."
    @assert length(action) == length(x.materials)*ne(x.network) "Reorder action vector must have length num_products * num_edges."
    act = show_action(x, action)
    #increase period counter
    x.period += 1
    #intialize next period on-hand and pipeline inventories with previous inventories
    initialize_inventories!(x)
    #move active shipments forward one period
    x.shipments.lead .-= 1
    #decrease time until due date for active orders by one period
    x.open_orders.due .-= 1
    #place requests
    place_orders!(x, act)
    #update on hand and pipeline inventories due to arrived shipments
    arrivals = update_shipments!(x)
    #discard any excess inventory
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)
    #markets open and demand occurs
    simulate_markets!(x)
    #update inventories at each node and echelon
    update_inventories!(x)
    #calculate profit at each node
    if x.options[:evaluate_profit]
        calculate_profit!(x, arrivals)
        x.reward = sum(filter(:period => j -> j == x.period, x.profit, view=true).value) #update reward (current profit). NOTE: is this ok for RL?
    end
end

"""
    initialize_inventories!(x::SupplyChainEnv)

Initialize on-hand and pipeline inventories with those from the previous period.
"""
function initialize_inventories!(x::SupplyChainEnv)
    #find previous period's inventories
    prev_inventory_on_hand = filter(:period => j -> j == x.period-1, x.inventory_on_hand, view=true) #previous inventory level
    prev_onhand_grp = groupby(prev_inventory_on_hand, [:node, :material])
    prev_inventory_pipeline = filter(:period => j -> j == x.period-1, x.inventory_pipeline, view=true) #previous inventory level
    prev_pipeln_grp = groupby(prev_inventory_pipeline, [:arc, :material])
    #update dataframes
    for mat in x.materials
        for n in vertices(x.network)
            prev = prev_onhand_grp[(node = n, material = mat)][1,:level] #previous on hand inventory
            push!(x.inventory_on_hand, [x.period, n, mat, prev, 0]) #intialize with previous inventory levels
        end
        for a in edges(x.network)
            prev = prev_pipeln_grp[(arc = (a.src,a.dst), material = mat)][1,:level] #previous pipeline inventory
            push!(x.inventory_pipeline, [x.period, (a.src,a.dst), mat, prev]) #intialize with previous inventory levels
        end
    end
end

"""
    place_orders!(x::SupplyChainEnv, act::NamedArray)

Place inventory replenishment requests throughout the network.
"""
function place_orders!(x::SupplyChainEnv, act::NamedArray)
    arcs = names(act,2) #arcs
    #exit if no action
    if iszero(act)
        exit_order!(x, arcs)
        return
    end
    #extract info
    mats = x.materials #materials
    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    #sample lead times
    leads = Dict((a,mat) => rand(get_prop(x.network, a, :lead_time)[mat]) for a in edges(x.network), mat in mats)
    #identify nodes that can place requests
    nodes = topological_sort_by_dfs(x.network) #sort nodes in topological order so that orders are placed moving down the network
    source_nodes = filter(n -> isempty(inneighbors(x.network, n)), nodes) #source nodes (can't place replenishment orders)
    request_nodes = setdiff(nodes, source_nodes) #nodes placing requests (all non-source nodes)
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inventory_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inventory_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #place requests
    for req in request_nodes #loop by nodes placing requests (in reverse topological order)
        sup_priority = get_prop(x.network, req, :supplier_priority) #get supplier priority list
        for mat in mats, sup in sup_priority[mat] #loop by material and supplier priority
            a = (sup, req) #arc
            lead = float(leads[Edge(a...), mat]) #sampled lead time
            amount = act[mat, a] #amount requested
            #calculate backlog (will be 0 if assuming lost sales)
            backlog = calculate_backlog(x, sup, req, mat)
            #continue to next iteration if no request made and no backlog
            if iszero(amount) && iszero(backlog)
                push!(x.demand, [x.period, a, mat, 0, 0, 0, 0, missing])
                continue 
            end
            #create order and save service time
            create_order!(x, a..., mat, amount)
            #try to fulfill order from stock
            fulfill_from_stock!(x, a..., mat, lead, supply_grp, pipeline_grp)
            #try to fulfill order from production
            if sup in x.producers
                fulfill_from_production!(x, a..., mat, lead, supply_grp, pipeline_grp, capacities)
            end
            #enforce guaranteed_service operation
            if x.options[:guaranteed_service] #any open orders with non-positive relative due date (0 service time) become lost sales
                filter!(:due => t -> t > 0, x.open_orders)
            end
        end
    end

    #updated production capacities (after commited production)
    for n in x.producers
        set_prop!(x.network, n, :production_capacity, capacities[n])
    end
end

"""
    exit_order!(x, arcs)

Abort order placement.
"""
function exit_order!(x::SupplyChainEnv, arcs::Vector)
    for a in arcs, mat in x.materials
        backlog = calculate_backlog(x, a..., mat)
        push!(x.demand, [x.period, a, mat, 0, 0, 0, backlog, missing])
    end
end

"""
    calculate_backlog(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol)

Calculate quantity of material `mat` backlogged on the arc (`sup`,`req`) (includes open orders).
"""
function calculate_backlog(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol)
    due_orders = filter([:arc, :material] => (a,m) -> a == (sup,req) && m == mat, x.open_orders, view=true)
    backlog = reduce(+, due_orders.quantity; init = 0)
    
    return backlog
end

"""
    create_order!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, amount::Real)

Create and log order.
"""
function create_order!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, amount::Real)
    x.num_orders += 1 #create new order ID
    push!(x.orders, [x.num_orders, x.period, (sup,req), mat, amount, []]) #update order history
    if sup == req #get service time for market node
        service_time = get_prop(x.network, sup, :service_time)[mat] 
    else #get service time for arc
        service_time = get_prop(x.network, sup, req, :service_time)[mat] 
    end
    push!(x.open_orders, [x.num_orders, (sup,req), mat, amount, service_time]) #add order to temp order df
    sort!(x.open_orders, [:due, :id]) #sort by service time and creation date (serve orders with lowest service time, ranked by order age)
end

"""
    fulfill_from_stock!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{Missing,GroupedDataFrame}
    )

Fulfill request from on-hand inventory.
"""
function fulfill_from_stock!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{Missing,GroupedDataFrame}
)
    #available supply
    supply = supply_grp[(node = src, material = mat)].level[1]
    
    #loop through orders
    orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.open_orders, view = true)
    for row in eachrow(orders_df)
        order_amount = row.quantity
        accepted_inv = min(order_amount, supply) #amount fulfilled from inventory
        if accepted_inv > 0
            row.quantity -= accepted_inv #update x.open_orders (deduct fulfilled part of the order)
            push!(x.orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_inv)) #update fulfilled column in order history (date, supplier, amount fulfilled)
            push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_inv, lead, 0, missing]) #log demand
            supply_grp[(node = src, material = mat)].level[1] -= accepted_inv #remove inventory from site
            !ismissing(pipeline_grp) && make_shipment!(x, src, dst, mat, accepted_inv, lead, pipeline_grp) #ship material (unless it is external demand)
        end
        if row.quantity > 0 && row.due <= 0 && !isproduced(x, src, mat) #if some amount of the order is due and wasn't fulfilled log it as unfulfilled & try to reallocate (only if node is not a plant that produces mat since the plant will try to fulfill from production next)
            log_unfulfilled_demand!(x, src, dst, mat, accepted_inv, row)
        end
    end
    #remove any fulfilled orders from x.open_orders
    filter!(:quantity => q -> q > 0, x.open_orders) 
end

"""
    fulfill_from_production!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
    )

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
)
    #extract info
    !isproduced(x, src, mat) && return #exit if mat is not produced at this plant
    bom = get_prop(x.network, src, :bill_of_materials)
    rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
    cmat_names = names(filter(k -> k > 0, bom[:,mat]), 1) #names of co-products

    #loop through orders
    orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.open_orders, view = true)
    for row in eachrow(orders_df)
        cap_and_sup = get_capacity_and_supply(src, mat, bom, rmat_names, cmat_names, capacities, supply_grp)
        order_amount = row.quantity #amount requested in order
        accepted_prod = min(order_amount, cap_and_sup...) #fulfill remaining with available capacity
        if accepted_prod > 0
            row.quantity -= accepted_prod #update x.open_orders (deduct fulfilled quantity)
            push!(x.orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_prod)) #update fulfilled column in order history (date, supplier, and amount fulfilled)
            push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_prod, lead, 0, missing]) #log demand
            capacities[src][mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
            for rmat in rmat_names #reactant consumed
                supply_grp[(node = src, material = rmat)].level[1] += accepted_prod * bom[rmat,mat]
            end
            make_shipment!(x, src, dst, mat, accepted_prod, lead, pipeline_grp) #schedule shipment
            for cmat in cmat_names #coproducts scheduled for production
                coproduction = accepted_prod*bom[cmat,mat]
                capacities[src][cmat] -= coproduction #update coproduct production capacity
                push!(x.demand, [x.period, (src,dst), cmat, 0, coproduction, lead, 0, missing]) #log demand
                make_shipment!(x, src, dst, cmat, coproduction, lead, pipeline_grp) #schedule shipment of co-product (original order amount is 0 since this is a coproduct)
            end
        end
        if row.quantity > 0 && row.due <= 0 #if some amount of the order is due and wasn't fulfilled, log it as unfulfilled & try to reallocate it
            log_unfulfilled_demand!(x, src, dst, mat, accepted_prod, row)
        end
    end
    #remove any fulfilled orders from x.open_orders
    filter!(:quantity => q -> q > 0, x.open_orders) 
end

"""
    isproduced(x::SupplyChainEnv, n::Int, mat::Symbol)

Check if material `mat` is produced in node `n`.
"""
function isproduced(x::SupplyChainEnv, n::Int, mat::Symbol)
    !in(:bill_of_materials, keys(props(x.network, n))) && return false #n is not a plant
    bom = get_prop(x.network, n, :bill_of_materials)
    !in(mat, names(bom,2)) && return false #mat is not produced at this plant
    raws = filter(k -> k < 0, bom[:,mat]) #names of raw materials
    isempty(raws) && return false #mat is not produced at this plant (no raw materials are converted to mat)

    return true
end

"""
    get_capacity_and_supply(
        bom::NamedArray, rmat_names::Vector, cmat_names::Vector, 
        capacities::Dict, supply_grp::GroupedDataFrame
    )

Get available capacity and material supply at producer.
"""
function get_capacity_and_supply(
    n::Int, mat::Symbol, bom::NamedArray, 
    rmat_names::Vector, cmat_names::Vector, 
    capacities::Dict, supply_grp::GroupedDataFrame
)
    #commit production at plant
    capacity = [] #get production capacity
    mat_supply = [] #store max capacity based on raw material consumption for each raw material
    isempty(rmat_names) && return [0],[0] #if material is not produced at the node, then return zero capacity and zero supply
    push!(capacity, capacities[n][mat]) #production capacity for that material
    for rmat in rmat_names #check raw material supply
        sup_pp = supply_grp[(node = n, material = rmat)].level[1] #supply of material involved in BOM
        push!(mat_supply, - sup_pp / bom[rmat,mat]) #only account for raw materials that are in the BOM
    end 
    for cmat in cmat_names #add capacity constraint for any co-products (scaled by stoichiometry)
        push!(capacity, capacities[n][rmat] / bom[cmat,mat])
    end

    return vcat(capacity, mat_supply)
end

"""
    make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)

Schedule shipment of material.
"""
function make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)
    #ship material
    push!(x.shipments, [(src,dst), mat, accepted_amount, lead])
    pipeline_grp[(arc = (src, dst), material = mat)].level[1] += accepted_amount #add inventory to pipeline
    if iszero(lead) #if zero leadtime, update shipments so that it is immediately available for fulfilling downstream orders
        update_shipments!(x)
    end
end

"""
    log_unfulfilled_demand!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, accepted::Float64, order_row::DataFrameRow)

Calculate unfulfilled demand and reallocate order to next supplier in the priority list.
"""
function log_unfulfilled_demand!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, accepted::Float64, order_row::DataFrameRow)
    new_alloc = missing
    if x.options[:reallocate] && sup != req #don't reallocate if external demand
        sup_priority = get_prop(x.network, req, :supplier_priority)[mat]
        sup_priority = vcat(sup_priority, sup_priority[1]) #add top priority to the end of the list so that lowest priority reallocates to top priority (cycle)
        next_sup_loc = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
        new_sup = sup_priority[next_sup_loc] #next supplier in line
        new_alloc = (new_sup, req) #store new arc where reallocated
        order_row.arc = new_alloc #update x.open_orders to reallocate material
    end
    if accepted > 0 #order was partially fulfilled
        x.demand[end,:unfulfilled] = order_row.quantity
        x.demand[end,:reallocated] = new_alloc
    else #order was not fulfilled
        original_amount = order_row.quantity
        push!(x.demand, [x.period, (sup,req), mat, original_amount, 0, 0, original_amount, new_alloc])
    end
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    #filter data
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inventory_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inventory_pipeline, view=true)
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
    on_hand_df = filter(:period => j -> j == x.period, x.inventory_on_hand, view=true) #on_hand inventories
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
    on_hand_df = filter(:period => j -> j == x.period, x.inventory_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    pipeline_df = filter(:period => j -> j == x.period, x.inventory_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, :material) #group by material
    orders_df = filter([:period, :reallocated] => (i1, i2) -> i1 == x.period && ismissing(i2), x.demand, view=true) #replenishment orders from current period
    orders_grp = groupby(orders_df, :material) #group by material

    #initialize echelon inventory positions
    initialize_echelons!(x)
    ech_df = filter(:period => j -> j == x.period, x.echelon_stock, view=true) #echelon stock at each node
    ech_grp = groupby(ech_df, [:node, :material]) #group by material
    
    #loop through nodes and update inventory levels, positions, and echelons
    for n in vertices(x.network), mat in x.materials
        ilevel, ipos = inventory_components(x, n, mat, pipeline_grp, onhand_grp, orders_grp)
        push!(x.inventory_level, [x.period, n, mat, ilevel]) #update inventory level
        push!(x.inventory_position, [x.period, n, mat, ipos]) #update inventory position
        update_echelons!(x, n, mat, ipos, ech_grp) #update echelon stocks
    end

end

"""
    initialize_echelons!(x::SupplyChainEnv)

Initialize echelon stocks at the current period to 0.
"""
function initialize_echelons!(x::SupplyChainEnv)
    for n in vertices(x.network), mat in x.materials
        if get_prop(x.network, n, :inventory_capacity)[mat] > 0
            push!(x.echelon_stock, [x.period, n, mat, 0])
        end
    end
end

"""
    inventory_components(
        x::SupplyChainEnv, n::Int, mat::Symbol, 
        pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, 
        orders_grp::GroupedDataFrame
    )

Extract components to determine inventory level and position.
"""
function inventory_components(
    x::SupplyChainEnv, n::Int, mat::Symbol, 
    pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, 
    orders_grp::GroupedDataFrame
)
    pipeline = sum(filter(:arc => j -> j[end] == n, pipeline_grp[(material = mat,)], view=true).level) #in-transit inventory
    onhand = onhand_grp[(node = n, material = mat)].level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        backlog = sum(filter(:arc => i -> i[1] == n, orders_grp[(material = mat,)], view=true).unfulfilled)
        backorder = sum(filter(:arc => i -> i[2] == n && i[2] != i[1], orders_grp[(material = mat,)], view=true).unfulfilled)
    end
    ilevel = onhand - backlog #inventory level
    iorder = pipeline + backorder #inventory on order
    ipos = ilevel + iorder #inventory position

    return ilevel, ipos
end

"""
    update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ech_grp::GroupedDataFrame)

Update echelon stocks for current time period.
"""
function update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ech_grp::GroupedDataFrame)
    for ech in findall(i -> n in i, x.echelons) #identify which echelons have been affected and add to these
        if get_prop(x.network, ech, :inventory_capacity)[mat] > 0 #only add to echelon if that node holds that material
            ech_grp[(node = ech, material = mat)].level[1] += ipos
        end
    end
end

"""
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    supply_df = filter([:period,:node] => (t,n) -> t == x.period && n in x.markets, x.inventory_on_hand, view=true) #on_hand inventory at node
    supply_grp = groupby(supply_df, [:node, :material]) #group by node and material
    for n in x.markets
        dmnd_seq_dict = get_prop(x.network, n, :demand_sequence)
        dmnd_freq_dict = get_prop(x.network, n, :demand_period)
        dmnd_dist_dict = get_prop(x.network, n, :demand_distribution)
        for mat in x.materials
            dmnd_seq = dmnd_seq_dict[mat]
            dprob = Bernoulli(1/dmnd_freq_dict[mat]) #demand probability (probability of ordering is 1/demand period; aka, once every x days)
            dmnd = dmnd_dist_dict[mat] #demand distribution
            q = iszero(dmnd_seq) ? rand(dprob) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
            create_order!(x, n, n, mat, q)
            fulfill_from_stock!(x, n, n, mat, 0., supply_grp, missing)
        end
    end
end