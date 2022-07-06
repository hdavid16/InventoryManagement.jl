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

    #check if partial fulfillment is allowed
    partial_fulfillment = check_fulfillment_type(x.network, src, dst, mat)
    
    #loop through orders on relevant arc
    orders_df = relevant_orders(x, src, dst, mat)
    for row in eachrow(orders_df)
        order_amount = row.quantity
        accepted_inv = accepted_from_stock(order_amount, supply, partial_fulfillment)
        if accepted_inv > 0
            row.quantity -= accepted_inv #update x.open_orders (deduct fulfilled part of the order)
            push!(x.fulfillments, (row.id, x.period, src, accepted_inv)) #update order fulfillments log (date, supplier, amount fulfilled)
            push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_inv, lead, 0, missing]) #log demand
            supply_grp[(node = src, material = mat)].level[1] -= accepted_inv #remove inventory from site
            dst != :market && make_shipment!(x, src, dst, mat, accepted_inv, lead, supply_grp, pipeline_grp) #ship material (unless it is external demand)
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
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{GroupedDataFrame,Missing}, capacities::Dict
    )

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(
    x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::Union{GroupedDataFrame,Missing}, capacities::Dict
)

    #check if partial fulfillment is allowed
    partial_fulfillment = check_fulfillment_type(x.network, src, dst, mat)

    #extract info
    bom = get_prop(x.network, src, :bill_of_materials)
    rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
    cmat_names = names(filter(k -> k > 0, bom[:,mat]), 1) #names of co-products

    #get relevant orders on that arc & raw material orders (:production)
    orders_df = relevant_orders(x, src, dst, mat) #NOTE: Assume that if plant requests production, it will accept early fulfillment
    raw_orders_grp = raw_material_orders(x, orders_df.id, rmat_names)
    #loop through orders on relevant arc
    for row in eachrow(orders_df)
        cap_and_sup = get_capacity_and_supply(src, mat, bom, rmat_names, cmat_names, capacities, supply_grp)
        order_amount = row.quantity #amount requested in order
        accepted_prod = accepted_production(order_amount, cap_and_sup, partial_fulfillment) #amount accepted
        if accepted_prod > 0
            #fulfill order (may be partial)
            row.quantity -= accepted_prod #update x.open_orders (deduct fulfilled quantity)
            #consume reactant
            for rmat in rmat_names
                consume_reactant!(row.id, src, rmat, bom[rmat,mat], accepted_prod, supply_grp, raw_orders_grp) 
            end
            #schedule coproduction
            for cmat in cmat_names
                co_production!(x, src, dst, cmat, bom[cmat,mat], accepted_prod, lead, supply_grp, pipeline_grp, capacities)
            end
            #schedule production
            production!(x, row.id, src, dst, mat, order_amount, accepted_prod, lead, supply_grp, pipeline_grp, capacities)
        end
        if row.quantity > 0 && row.due <= 0 #if some amount of the order is due and wasn't fulfilled, log it as unfulfilled & try to reallocate it
            log_unfulfilled_demand!(x, row, accepted_prod)
        end
    end
    #remove any fulfilled orders from x.open_orders
    filter!(:quantity => q -> q > 0, x.open_orders) 
end

"""
    check_fulfillment_type(net::MetaDiGraph, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})

Check if destination `dst` node accepts partial fulfillments.
"""
function check_fulfillment_type(net::MetaDiGraph, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})
    #check if partial fulfillment is allowed
    src == dst && return true #NOTE: ASSUME if self-request (i.e., production order), then allow partial fulfillments
    indicator_node = dst == :market ? src : dst
    param_key = dst == :market ? :market_partial_fulfillment : :partial_fulfillment
    partial_fulfillment = get_prop(net, indicator_node, param_key)[mat]

    return partial_fulfillment
end

"""
    accepted_from_stock(order_amount, supply, partial_fulfillment)

Amount accepted from available supply on an order received.
"""
function accepted_from_stock(order_amount, supply, partial_fulfillment)
    if partial_fulfillment
        return min(order_amount, supply) #amount fulfilled from inventory
    else
        return order_amount <= supply ? order_amount : 0. #only accept full order or nothing at all
    end
end

"""
    accepted_production(order_amount, cap_and_sup, partial_fulfillment)

Amount of production order accepted based on available capacity, supply, and partial fulfillment criteria.
"""
function accepted_production(order_amount, cap_and_sup, partial_fulfillment)
    if partial_fulfillment
        return min(order_amount, cap_and_sup...) #fulfill remaining with available capacity
    else
        return order_amount <= minimum(cap_and_sup) ? order_amount : 0. #only accept full order or nothing at all
    end
end

"""
    consume_reactant!(order_id, src, rmat, stoich, accepted_prod, supply_grp, raw_orders_grp)

Consume reactant `rmat` and update inventory and orders.
"""
function consume_reactant!(order_id, src, rmat, stoich, accepted_prod, supply_grp, raw_orders_grp)
    consumed = accepted_prod * stoich #negative number
    supply_grp[(node = src, material = rmat)].level[1] += consumed
    raw_orders_grp[(id = order_id, material = rmat)].quantity[1] += consumed
end

"""
    co_production!(x::SupplyChainEnv, src, dst, cmat, stoich, accepted_prod, lead, supply_grp, pipeline_grp, capacities)

Produce coproduct `cmat` and update inventory and shipments
"""
function co_production!(x::SupplyChainEnv, src, dst, cmat, stoich, accepted_prod, lead, supply_grp, pipeline_grp, capacities)
    coproduction = accepted_prod*stoich
    capacities[cmat] -= coproduction #update coproduct production capacity
    push!(x.demand, [x.period, (src,dst), cmat, 0, coproduction, lead, 0, missing]) #log demand
    if dst == :market && iszero(lead) #only triggered if make-to-order & has 0 production lead time
        supply_grp[(node = src, material = cmat)].level[1] += coproduction
    else
        make_shipment!(x, src, dst, cmat, coproduction, lead, supply_grp, pipeline_grp) #schedule shipment of co-product (original order amount is 0 since this is a coproduct)
    end
end

"""
    production!(x::SupplyChainEnv, src, dst, mat, order_amount, accepted_prod, lead, supply_grp, pipeline_grp, capacities)

Produce material `mat` and create shipments and update fulfillments/demand tables.
"""
function production!(x::SupplyChainEnv, order_id, src, dst, mat, order_amount, accepted_prod, lead, supply_grp, pipeline_grp, capacities)
    capacities[mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
    push!(x.fulfillments, (order_id, x.period, src, accepted_prod)) #update order fulfillments (date, supplier, and amount fulfilled)
    push!(x.demand, [x.period, (src,dst), mat, order_amount, accepted_prod, lead, 0, missing]) #log demand
    dst != :market && make_shipment!(x, src, dst, mat, accepted_prod, lead, supply_grp, pipeline_grp) 
end

"""
    make_shipment!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, 
        accepted_amount::Float64, lead::Float64, 
        supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame
    )

Schedule shipment of material.
"""
function make_shipment!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Union{Symbol,String}, 
    accepted_amount::Float64, lead::Float64, 
    supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame
)
    #add inventory to pipeline
    a = (src, dst)
    pipeline_grp[(arc = a, material = mat)].level[1] += accepted_amount 
    if iszero(lead) #if zero leadtime, update shipments so that it is immediately available for fulfilling downstream orders
        shipment_completed!(x, accepted_amount, a, mat, supply_grp, pipeline_grp)
    else #otherwise, log the shipment
        push!(x.shipments, [a, mat, accepted_amount, lead])
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
        #update inventories and capacities
        shipment_completed!(x, amount, a, mat, supply_grp, pipeline_grp)
    end
    filter!(:lead => i -> i > 0, x.shipments) #remove shipments that arrived (and any zero values)
    
    return arrivals
end

"""
    shipment_completed!(x::SupplyChainEnv, amount::Float64, arc::Tuple, mat::Union{Symbol,String}, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)

Update downstream inventory and pipeline inventory for completed shipment, and update production capacities for completed production.
"""
function shipment_completed!(x::SupplyChainEnv, amount::Float64, arc::Tuple, mat::Union{Symbol,String}, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame)
    #update downstream inventory and pipeline inventory
    supply_grp[(node = arc[end], material = mat)].level[1] += amount
    pipeline_grp[(arc = arc, material = mat)].level[1] -= amount
    #restore production capacities (production has concluded)
    if arc[1] == arc[2]
        #restore mat production capacity
        capacity = get_prop(x.network, arc[1], :production_capacity)
        capacity[mat] += amount
        #restore any co-product production capacity
        bom = get_prop(x.network, arc[1], :bill_of_materials)
        cmat_names = names(filter(k -> k > 0, bom[:,mat]), 1) #names of co-products
        for cmat in cmat_names
            capacity[cmat] += amount*bom[cmat,mat]
        end
        #update production capacities
        set_prop!(x.network, arc[1], :production_capacity, capacity)
    end
end