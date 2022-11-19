"""
    fulfill_from_stock!(
        x::SupplyChainEnv, src::Int, dst::Vector{Int}, mat::Material, leads::Dict
    )

Fulfill request from on-hand inventory to downstream nodes.

    fulfill_from_stock!(x::SupplyChainEnv, src::Int, dst::Symbol, mat::Material)

Fulfill market (`dst`) request from on-hand inventory at `src`.
"""
function fulfill_from_stock!(
    x::SupplyChainEnv, src::Int, dsts::Vector{Int}, mat::Material, leads::Dict
)    
    #loop through orders on relevant arc
    orders_df = relevant_orders(x, src, dsts, mat)
    reallocated_suppliers = [] #store reallocated orders
    for row in eachrow(orders_df)
        supply = x.tmp[src,mat,:on_hand] #available supply
        partial_fulfillment = get_prop(x.network, row.arc[2], :partial_fulfillment)[mat] #check if partial fulfillment is allowed
        accepted_inv = accepted_from_stock(row.amount, supply, partial_fulfillment)
        if accepted_inv > 0
            #log fulfillment
            row.amount = floor(row.amount - accepted_inv, digits=x.options[:numerical_precision]) #update x.open_orders (deduct fulfilled part of the order)
            x.tmp[src,mat,:on_hand] -= accepted_inv #remove inventory from site
            push!(x.fulfillments, (row.id, x.period, row.arc, mat, accepted_inv, :sent)) #log order as sent
            #ship material
            lead = leads[Edge(row.arc...), mat] #sampled lead time
            make_shipment!(x, row.id, row.arc..., mat, accepted_inv, lead) #ship material (unless it is external demand)
        end
        if x.options[:reallocate] #try to reallocate unfulfilled order if reallocation allowed
            new_sup = reallocate_demand!(x, row) #NOTE: New fulfillment logic doesn't allow reallocation to work properly
            if !isnothing(new_sup)
                push!(reallocated_suppliers, new_sup)
            end
        end
    end
    update_orders!(x, orders_df)

    return unique(reallocated_suppliers)
end
function fulfill_from_stock!(x::SupplyChainEnv, src::Int, dst::Symbol, mat::Material)
    #check if partial fulfillment is allowed
    partial_fulfillment = get_prop(x.network, src, :market_partial_fulfillment)[mat]
    
    #loop through orders on relevant arc
    orders_df = relevant_orders(x, src, dst, mat)
    for row in eachrow(orders_df)
        supply = x.tmp[src,mat,:on_hand] #available supply
        accepted_inv = accepted_from_stock(row.amount, supply, partial_fulfillment)
        if accepted_inv > 0
            #log fulfillment
            row.amount = floor(row.amount - accepted_inv, digits=x.options[:numerical_precision]) #update x.open_orders (deduct fulfilled part of the order)
            x.tmp[src,mat,:on_hand] -= accepted_inv #remove inventory from site
            push!(x.fulfillments, (row.id, x.period, row.arc, mat, accepted_inv, :sent)) #log order as sent
            push!(x.fulfillments, (row.id, x.period, row.arc, mat, accepted_inv, :delivered)) #log order as delivered
        end
    end
    update_orders!(x, orders_df)
end

"""
    fulfill_from_production!(
        x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material, 
        lead::Real, capacities::Dict
    )

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(
    x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material, 
    lead::Int, capacities::Dict
)
    #check if partial fulfillment is allowed
    partial_fulfillment = src == dst || get_prop(x.network, src, :market_partial_fulfillment) #NOTE: ASSUME if self-request (i.e., production order), then allow partial fulfillments

    #extract info
    bom = get_prop(x.network, src, :bill_of_materials)
    rmat_names = names(filter(<(0), bom[:,mat]), 1) #names of raw materials
    cmat_names = names(filter(>(0), bom[:,mat]), 1) #names of co-products

    #get relevant orders on that arc & raw material orders (:consumption)
    orders_df = relevant_orders(x, src, dst, mat) #NOTE: Assume that if plant requests production, it will accept early fulfillment
    raw_orders_grp = raw_material_orders(x, orders_df.id, rmat_names)
    #loop through orders on relevant arc
    for row in eachrow(orders_df)
        cap_and_sup = get_capacity_and_supply(x, src, mat, bom, rmat_names, cmat_names, capacities)
        accepted_prod = accepted_production(row.amount, cap_and_sup, partial_fulfillment) #amount accepted
        if accepted_prod > 0
            #fulfill order (may be partial)
            row.amount = floor(row.amount - accepted_prod, digits=x.options[:numerical_precision]) #update x.open_orders (deduct fulfilled quantity)
            #consume reactant
            for rmat in rmat_names
                consume_reactant!(x, row.id, src, rmat, bom[rmat,mat], accepted_prod, raw_orders_grp) 
            end
            #schedule coproduction
            for cmat in cmat_names
                co_production!(x, row.id, src, cmat, bom[cmat,mat], accepted_prod, lead, capacities)
            end
            #schedule production
            production!(x, row.id, src, dst, mat, accepted_prod, lead, capacities)
        end
    end
    #remove any fulfilled orders from x.open_orders
    fulfilled_orders = Set(filter(:amount => <=(0), orders_df, view=true).id) #find fulfilled orders
    if !isempty(fulfilled_orders)
        filter(:id => in(fulfilled_orders), x.orders, view=true).fulfilled .= x.period #store fulfillment date
        filter!(:id => !in(fulfilled_orders), x.open_orders) #remove fulfilled orders (with any associated production orders)
    end
end

# """
#     check_fulfillment_type(net::MetaDiGraph, src::Int, dst::Union{Int,Symbol}, mat::Material)

# Check if destination `dst` node accepts partial fulfillments.
# """
# check_fulfillment_type(net::MetaDiGraph, src::Int, dst::Int, mat::Material) = get_prop(net, dst, :partial_fulfillment)[mat]
# check_fulfillment_type(net::MetaDiGraph, src::Int, dst::Symbol, mat::Material) = get_prop(net, src, :market_partial_fulfillment)[mat]

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
    consume_reactant!(x::SupplyChainEnv, order_id, src, rmat, stoich, accepted_prod, raw_orders_grp)

Consume reactant `rmat` and update inventory and orders.
"""
function consume_reactant!(x::SupplyChainEnv, order_id, src, rmat, stoich, accepted_prod, raw_orders_grp)
    consumed = accepted_prod * stoich #negative number
    x.tmp[src,rmat,:on_hand] += consumed
    raw_orders_grp[(id = order_id, material = rmat)].amount[1] += consumed
    push!(x.fulfillments, (order_id, x.period, (src,:consumption), rmat, -consumed, :sent)) #raw material is consumed ("sent")
end

"""
    co_production!(x::SupplyChainEnv, order_id, src, cmat, stoich, accepted_prod, lead, capacities)

Produce coproduct `cmat` and update inventory and shipments
"""
function co_production!(x::SupplyChainEnv, order_id, src, cmat, stoich, accepted_prod, lead, capacities)
    coproduction = accepted_prod*stoich
    capacities[cmat] -= coproduction #update coproduct production capacity
    push!(x.fulfillments, (order_id, x.period, (src,:coproduction), cmat, coproduction, :sent))
    make_shipment!(x, order_id, src, src, cmat, coproduction, lead) #schedule shipment of co-product 
end

"""
    production!(x::SupplyChainEnv, src, dst, mat, accepted_prod, lead, capacities)

Produce material `mat` and create shipments and update fulfillments/demand tables.
"""
function production!(x::SupplyChainEnv, order_id, src, dst, mat, accepted_prod, lead, capacities)
    capacities[mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
    push!(x.fulfillments, (order_id, x.period, (src,dst), mat, accepted_prod, :sent)) #update fulfillment log
    make_shipment!(x, order_id, src, dst, mat, accepted_prod, lead) 
end

"""
    make_shipment!(
        x::SupplyChainEnv, order_id::Int,
        src::Int, dst::Union{Int,Symbol}, mat::Material, 
        accepted_amount::Float64, lead::Float64
    )

Schedule shipment of material.
"""
function make_shipment!(
    x::SupplyChainEnv, order_id::Int,
    src::Int, dst::Union{Int,Symbol}, mat::Material, 
    accepted_amount::Float64, lead::Int
)
    #add inventory to pipeline
    a = (src,dst)
    if dst != :market #don't add if market order
        x.tmp[a,mat,:pipeline] += accepted_amount
    end
    push!(x.shipments, [order_id, a, mat, accepted_amount, lead]) #log the shipment
    iszero(lead) && update_shipments!(x) #if leadtime is zero, update shipments immediately
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    #find active shipments with 0 lead time
    arrivals = filter(:lead => <=(0), x.shipments, view=true) 
    for row in eachrow(arrivals)
        id, a, mat, amount = row[[:id, :arc, :material, :amount]]
        #update inventories and capacities
        shipment_completed!(x, id, amount, a, mat)
    end
    filter!(:lead => >(0), x.shipments) #remove shipments that arrived (and any zero values)
end

"""
    shipment_completed!(x::SupplyChainEnv, order_id::Int, amount::Float64, arc::Tuple, mat::Material)

Update downstream inventory and pipeline inventory for completed shipment, and update production capacities for completed production.
"""
function shipment_completed!(x::SupplyChainEnv, order_id::Int, amount::Float64, arc::Tuple, mat::Material)
    #update downstream inventory and pipeline inventory
    if arc[end] != :market
        x.tmp[arc[end],mat,:on_hand] += amount
        x.tmp[arc,mat,:pipeline] -= amount
    end
    push!(x.fulfillments, (order_id, x.period, arc, mat, amount, :delivered))
    #restore production capacities (production has concluded)
    if arc[1] == arc[2] || arc[2] == :market
        #restore mat production capacity
        capacity = get_prop(x.network, arc[1], :production_capacity)
        capacity[mat] += amount
        #restore any co-product production capacity
        bom = get_prop(x.network, arc[1], :bill_of_materials)
        cmat_names = names(filter(>(0), bom[:,mat]), 1) #names of co-products
        for cmat in cmat_names
            capacity[cmat] += amount*bom[cmat,mat]
        end
        #update production capacities
        set_prop!(x.network, arc[1], :production_capacity, capacity)
    end
end