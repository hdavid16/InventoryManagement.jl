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
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    on_hand_df = filter(:period => j -> j == x.period, x.inventory_on_hand, view=true) #on_hand inventories
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material for easier lookup
    for n in vertices(x.network)
        node_mats = get_prop(x.network, n, :node_materials)
        node_cap = get_prop(x.network, n, :inventory_capacity)
        for mat in node_mats
            max_inv = node_cap[mat]
            isinf(max_inv) && continue #skip iteration if node doesn't store that material or has infinite capacity
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
    pipeline_grp = groupby(pipeline_df, [:arc, :material]) #group by arc and material
    due_by = x.options[:adjusted_stock] ? Inf : 0 #Inf means that all orders placed are counted (even if not due); otherwise, only due orders are counted
    orders_df = filter(:due => j -> j <= due_by, x.open_orders, view=true) #orders to count in backlogging
    orders_grp = groupby(orders_df, [:arc, :material]) #group by arc and material

    #initialize echelon inventory positions
    initialize_echelons!(x)
    ech_df = filter(:period => j -> j == x.period, x.echelon_stock, view=true) #echelon stock at each node
    ech_grp = groupby(ech_df, [:node, :material]) #group by material
    
    #loop through nodes and update inventory levels, positions, and echelons
    for n in vertices(x.network)
        for mat in get_prop(x.network, n, :node_materials) #loop through materials that can be stored at that node #x.materials
            ilevel, ipos = inventory_components(x, n, mat, pipeline_grp, onhand_grp, orders_grp)
            push!(x.inventory_level, [x.period, n, mat, ilevel]) #update inventory level
            push!(x.inventory_position, [x.period, n, mat, ipos]) #update inventory position
            update_echelons!(x, n, mat, ipos, ech_grp) #update echelon stocks
        end
    end

end

"""
    calculate_pending(arcs::Vector, mat::Union{Symbol,String}, orders_grp::GroupedDataFrame)

Calculate commited material quantity (open orders) using a pre-filtered grouped dataframe.
"""
calculate_pending(arcs::Vector, mat::Union{Symbol,String}, orders_grp::GroupedDataFrame) =
    sum(
        sum(orders_grp[(arc = a, material = mat)].quantity; init=0) 
        for a in arcs if (arc = a, material = mat) in keys(orders_grp); 
        init = 0
    )

"""
    calculate_in_transit(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, pipeline_grp::GroupedDataFrame)

Calculate in-transit inventory of material `mat` to node `n`.
"""
calculate_in_transit(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, pipeline_grp::GroupedDataFrame) = 
    sum(
        sum(pipeline_grp[(arc = (src,n), material = mat)].level; init=0) 
        for src in inneighbors(x.network, n) if (arc = (src,n), material = mat) in keys(pipeline_grp); 
        init=0
    )

"""
    inventory_components(
        x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, 
        pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, 
        orders_grp::GroupedDataFrame
    )

Extract components to determine inventory level and position.
"""
function inventory_components(
    x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, 
    pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, orders_grp::GroupedDataFrame
)
    pipeline = calculate_in_transit(x,n,mat,pipeline_grp) #in-transit inventory to node n
    onhand = onhand_grp[(node = n, material = mat)].level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        arcs_out = vcat( #nodes accounted for in backlog
            (n,:production), #raw material conversion
            (n,:market), #market sales
            [(n,succ) for succ in outneighbors(x.network,n) if n != succ] #downstream replenishments
        )
        backlog = calculate_pending(arcs_out,mat,orders_grp) 
        arcs_in = [(pred,n) for pred in inneighbors(x.network, n)] #backorder includes previous replenishment orders placed to upstream nodes
        backorder = calculate_pending(arcs_in,mat,orders_grp) 
    end
    ilevel = onhand - backlog #inventory level
    iorder = pipeline + backorder #inventory on order
    ipos = ilevel + iorder #inventory position

    return ilevel, ipos
end


"""
    initialize_echelons!(x::SupplyChainEnv)

Initialize echelon stocks at the current period to 0.
"""
function initialize_echelons!(x::SupplyChainEnv)
    for n in vertices(x.network), mat in get_prop(x.network, n, :node_materials)
        push!(x.echelon_stock, [x.period, n, mat, 0])
    end
end

"""
    update_echelons!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, ipos::Float64, ech_grp::GroupedDataFrame)

Update echelon stocks for current time period.
"""
function update_echelons!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, ipos::Float64, ech_grp::GroupedDataFrame)
    for ech in findall(i -> n in i, x.echelons) #identify which echelons have been affected and add to these
        if mat in get_prop(x.network, ech, :node_materials) #only add to echelon if that node holds that material
            ech_grp[(node = ech, material = mat)].level[1] += ipos
        end
    end
end