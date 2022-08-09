"""
    initialize_inventories!(x::SupplyChainEnv)

Initialize on-hand and pipeline inventories with those from the previous period.
"""
function initialize_inventories!(x::SupplyChainEnv)
    #update current inventory levels
    for n in vertices(x.network)
        for mat in get_prop(x.network, n, :node_materials)
            x.tmp[n,mat,:echelon] = 0 #initialize current echelon stock with 0
            x.tmp[n,mat,:discarded] = 0 #no inventory discarded initially
        end
    end
end

"""
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    for n in vertices(x.network)
        node_mats = get_prop(x.network, n, :node_materials)
        node_cap = get_prop(x.network, n, :inventory_capacity)
        for mat in node_mats
            max_inv = node_cap[mat]
            isinf(max_inv) && continue #skip iteration if node doesn't store that material or has infinite capacity
            onhand = x.tmp[n,mat,:on_hand] #onhand inventory
            if onhand > max_inv
                x.tmp[n,mat,:on_hand] = max_inv
                x.tmp[n,mat,:discarded] += onhand - max_inv
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
    due_by = x.options[:adjusted_stock] ? Inf : 0 #Inf means that all orders placed are counted (even if not due); otherwise, only due orders are counted
    orders_grp = groupby(
        filter(:due => <=(due_by), x.open_orders, view=true), #orders to count in backlogging 
        [:arc, :material]
    )
    
    #loop through nodes and update inventory levels, positions, and echelons
    for n in vertices(x.network)
        for mat in get_prop(x.network, n, :node_materials) #loop through materials that can be stored at that node #x.materials
            inventory_components(x, n, mat, orders_grp)
            update_echelons!(x, n, mat) #update echelon stocks
        end
    end

    #update dataframes in x (on_hand, pipeline, echelon)
    update_dfs!(x)
end

"""
    calculate_pending(arcs::Vector, mat::Union{Symbol,String}, orders_grp::GroupedDataFrame)

Calculate commited material quantity (open orders) using a pre-filtered grouped dataframe.
"""
calculate_pending(arcs::Vector, mat::Union{Symbol,String}, orders_grp::GroupedDataFrame) =
    sum(
        sum(orders_grp[(arc = a, material = mat)].amount; init=0) 
        for a in arcs if (arc = a, material = mat) in keys(orders_grp); 
        init = 0
    )

"""
    calculate_in_transit(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})

Calculate in-transit inventory of material `mat` to node `n`.
"""
calculate_in_transit(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}) = 
    sum(
        x.tmp[(src,n),mat,:pipeline] for src in inneighbors(x.network, n) 
            if mat in get_prop(x.network, src, :node_materials);
        init=0
    )

"""
    inventory_components(
        x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, 
        orders_grp::GroupedDataFrame
    )

Extract components to determine inventory level and position.
"""
function inventory_components(
    x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, 
    orders_grp::GroupedDataFrame
)
    pipeline = calculate_in_transit(x,n,mat) #in-transit inventory to node n
    onhand = x.tmp[n,mat,:on_hand] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        arcs_out = vcat( #nodes accounted for in backlog
            (n,:consumption), #raw material conversion
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

    x.tmp[n,mat,:level] = ilevel
    x.tmp[n,mat,:position] = ipos
end

"""
    update_echelons!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})

Update echelon stocks for current time period.
"""
function update_echelons!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})
    for ech in findall(i -> n in i, x.echelons) #identify which echelons have been affected and add to these
        if mat in get_prop(x.network, ech, :node_materials) #only add to echelon if that node holds that material
            x.tmp[ech,mat,:echelon] += x.tmp[n,mat,:position]
        end
    end
end

function update_dfs!(x::SupplyChainEnv)
    for n in vertices(x.network)
        for m in get_prop(x.network, n, :node_materials)
            #add rows to dfs
            if !iszero(x.tmp[n,m,:discarded])
                push!(x.inventory, (x.period, n, m, x.tmp[n,m,:discarded], :discarded))
            end
            push!(x.inventory, (x.period, n, m, x.tmp[n,m,:on_hand], :on_hand))
            push!(x.inventory, (x.period, n, m, x.tmp[n,m,:level], :level))
            push!(x.inventory, (x.period, n, m, x.tmp[n,m,:position], :position))
            push!(x.inventory, (x.period, n, m, x.tmp[n,m,:echelon], :echelon))
            for dst in outneighbors(x.network, n)
                arc = (n,dst)
                push!(x.inventory, (x.period, arc, m, x.tmp[arc,m,:pipeline], :pipeline))
            end
        end
    end
end