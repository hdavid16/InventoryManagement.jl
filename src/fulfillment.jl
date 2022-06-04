"""
    fulfill_from_stock!(
        x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{Missing,GroupedDataFrame}
    )

Fulfill request from on-hand inventory.
"""
function fulfill_from_stock!(
    x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{Missing,GroupedDataFrame}
)
    #available supply
    supply = supply_grp[(node = src, material = mat)].level[1]
    #node from which to check early_fulfillment and partial_fulfillment param values
    indicator_node = dst == :market ? src : dst
    
    #loop through orders
    if get_prop(x.network, indicator_node, :early_fulfillment)[mat]
        orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.open_orders, view = true)
    else #only attempt to fulfill orders that have expired their service lead time
        orders_df = filter([:arc, :material, :due] => (a,m,t) -> a == (src,dst) && m == mat && t <= 0, x.open_orders, view = true)
    end
    for row in eachrow(orders_df)
        order_amount = row.quantity
        if get_prop(x.network, indicator_node, :partial_fulfillment)[mat]
            accepted_inv = min(order_amount, supply) #amount fulfilled from inventory
        else
            accepted_inv = order_amount <= supply ? order_amount : 0 #only accept full order or nothing at all
        end
        if accepted_inv > 0
            row.quantity -= accepted_inv #update x.open_orders (deduct fulfilled part of the order)
            push!(x.orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_inv)) #update fulfilled column in order history (date, supplier, amount fulfilled)
            push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_inv, lead, 0, missing]) #log demand
            supply_grp[(node = src, material = mat)].level[1] -= accepted_inv #remove inventory from site
            dst != :market && make_shipment!(x, src, dst, mat, accepted_inv, lead, pipeline_grp) #ship material (unless it is external demand)
        end
        if row.quantity > 0 && row.due <= 0 #if some amount of the order is due and wasn't fulfilled log it as unfulfilled & try to reallocate
            log_unfulfilled_demand!(x, row, accepted_inv)
        end
    end
    #remove any fulfilled orders from x.open_orders
    filter!(:quantity => q -> q > 0, x.open_orders) 
end

"""
    fulfill_from_production!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
    )

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
)
    #extract info
    bom = get_prop(x.network, src, :bill_of_materials)
    rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
    cmat_names = names(filter(k -> k > 0, bom[:,mat]), 1) #names of co-products

    #loop through orders
    if get_prop(x.network, dst, :early_fulfillment)[mat]
        orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.open_orders, view = true)
    else #only attempt to fulfill orders that have expired their service lead time
        orders_df = filter([:arc, :material, :due] => (a,m,t) -> a == (src,dst) && m == mat && t <= 0, x.open_orders, view = true)
    end
    raw_orders_df = filter([:id, :material] => (id,m) -> id in orders_df.id && m != mat, x.open_orders, view = true)
    raw_orders_grp = groupby(raw_orders_df, [:id, :material])
    for row in eachrow(orders_df)
        cap_and_sup = get_capacity_and_supply(src, mat, bom, rmat_names, cmat_names, capacities, supply_grp)
        order_amount = row.quantity #amount requested in order
        if get_prop(x.network, dst, :partial_fulfillment)[mat]
            accepted_prod = min(order_amount, cap_and_sup...) #fulfill remaining with available capacity
        else
            accepted_prod = order_amount <= minimum(cap_and_sup) ? order_amount : 0 #only accept full order or nothing at all
        end
        if accepted_prod > 0
            #fulfill order (may be partial)
            row.quantity -= accepted_prod #update x.open_orders (deduct fulfilled quantity)
            #consume reactant
            for rmat in rmat_names 
                consumed = accepted_prod * bom[rmat,mat] #negative number
                supply_grp[(node = src, material = rmat)].level[1] += consumed
                raw_orders_grp[(id = row.id, material = rmat)].quantity[1] += consumed
            end
            #schedule coproduction
            for cmat in cmat_names
                coproduction = accepted_prod*bom[cmat,mat]
                capacities[src][cmat] -= coproduction #update coproduct production capacity
                push!(x.demand, [x.period, (src,dst), cmat, 0, coproduction, lead, 0, missing]) #log demand
                make_shipment!(x, src, dst, cmat, coproduction, lead, pipeline_grp) #schedule shipment of co-product (original order amount is 0 since this is a coproduct)
            end
            #schedule production
            capacities[src][mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
            push!(x.orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_prod)) #update fulfilled column in order history (date, supplier, and amount fulfilled)
            push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_prod, lead, 0, missing]) #log demand
            make_shipment!(x, src, dst, mat, accepted_prod, lead, pipeline_grp) 
        end
        if row.quantity > 0 && row.due <= 0 #if some amount of the order is due and wasn't fulfilled, log it as unfulfilled & try to reallocate it
            log_unfulfilled_demand!(x, row, accepted_prod)
        end
    end
    #remove any fulfilled orders from x.open_orders
    filter!(:quantity => q -> q > 0, x.open_orders) 
end

"""
    make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, accepted_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)

Schedule shipment of material.
"""
function make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, accepted_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)
    #ship material
    push!(x.shipments, [(src,dst), mat, accepted_amount, lead])
    pipeline_grp[(arc = (src, dst), material = mat)].level[1] += accepted_amount #add inventory to pipeline
    if iszero(lead) #if zero leadtime, update shipments so that it is immediately available for fulfilling downstream orders
        update_shipments!(x)
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