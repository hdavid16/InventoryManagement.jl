"""
    place_orders!(x::SupplyChainEnv, act::NamedArray)

Place inventory replenishment requests throughout the network.
"""
function place_orders!(x::SupplyChainEnv, act::NamedArray)
    #if no replenishments requested and no open orders exist, exit
    iszero(act) && isempty(x.open_orders) && return
    #extract info
    mats = x.materials #materials
    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    #sample lead times and service times
    leads = Dict((a,mat) => ceil(Int,rand(get_prop(x.network, a, :lead_time)[mat])) for a in edges(x.network), mat in mats)
    servs = Dict((a,mat) => ceil(Int,rand(get_prop(x.network, a, :service_lead_time)[mat])) for a in edges(x.network), mat in mats)
    #identify nodes that can place requests
    nodes = topological_sort(x.network) #sort nodes in topological order so that orders are placed moving down the network
    source_nodes = filter(n -> isempty(inneighbors(x.network, n)), nodes) #source nodes (can't place replenishment orders)
    request_nodes = setdiff(nodes, source_nodes) #nodes placing requests (all non-source nodes)
    #get copy of existing open orders
    orders_grp = groupby(copy(x.open_orders), [:arc, :material])

    #place requests
    for req in request_nodes #loop by nodes placing requests (in reverse topological order)
        sup_priority = get_prop(x.network, req, :supplier_priority) #get supplier priority list
        req_mats = get_prop(x.network, req, :node_materials) #get materials that are compatible with requester node
        for mat in req_mats, sup in sup_priority[mat] #loop by material and supplier priority
            a = (sup, req) #arc
            lead = leads[Edge(a...), mat] #sampled lead time
            serv = servs[Edge(a...), mat] #sampled service lead time
            amount = act[mat, a] #amount requested
            #continue to next iteration if no request made and no backlog
            if iszero(amount)
                #determine if there are any open orders for this material on the arc. If not, skip the fulfillment step.
                pending = calculate_pending([a,(req,:consumption)], mat, orders_grp) #include any raw material conversion orders (:consumption)
                if iszero(pending)
                    continue 
                end
            else
                #create order and save service lead time
                create_order!(x, a..., mat, amount, serv)
            end
            #try to fulfill request
            if sup == req
                fulfill_from_production!(x, a..., mat, lead, capacities[sup])
            else
                fulfill_from_stock!(x, a..., mat, lead)
            end
        end
    end

    #lost sales for expired orders (if backlog = false, or guaranteed service = true)
    if !x.options[:backlog] || x.options[:guaranteed_service] #any open orders with non-positive relative due date (0 service lead time) become lost sales
        lost_sales!(x)
    end

    #updated production capacities (after commited production)
    for n in x.producers
        set_prop!(x.network, n, :production_capacity, capacities[n])
    end
end

"""
    lost_sales!(x::SupplyChainEnv)

Delete expired orders (lost sales).
"""
function lost_sales!(x::SupplyChainEnv)
    expired_orders = @rsubset(x.open_orders, :due <= 0, view=true)
    for row in eachrow(expired_orders)
        push!(x.fulfillments, (row.id, x.period, row.arc, row.material, row.amount, :lost_sale))
    end
    if !isempty(expired_orders) #remove expired orders if there are any
        @rsubset!(x.open_orders, :due > 0)
    end
end

"""
    create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Float64)

Create and log order.
"""
function create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Real)
    x.num_orders += 1 #create new order ID
    push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,req), mat, amount]) #update order history
    push!(x.open_orders, [x.num_orders, service_lead_time, (sup,req), mat, amount]) #add order to temp order df
    #create raw material consumption and coproduction orders
    if (sup == req && isproduced(x,sup,mat)) || (req == :market && ismto(x,sup,mat))
        bom = get_prop(x.network, sup, :bill_of_materials)
        rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
        cmat_names = names(filter(k -> k > 0, bom[:,mat]), 1) #names of raw materials
        for rmat in rmat_names
            push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,:consumption), rmat, -amount*bom[rmat,mat]]) #update order history
            push!(x.open_orders, [x.num_orders, service_lead_time, (sup,:consumption), rmat, -amount*bom[rmat,mat]])
        end
        for cmat in cmat_names
            push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,:coproduction), cmat, amount*bom[cmat,mat]]) #update order history
            push!(x.open_orders, [x.num_orders, service_lead_time, (sup,:coproduction), cmat, amount*bom[cmat,mat]])
        end
    end
end

"""
    relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat`.
"""
function relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})
    #node from which to check early_fulfillment param value
    indicator_node = dst == :market ? src : dst
    param_key = dst == :market ? :market_early_fulfillment : :early_fulfillment
    early_fulfillment = get_prop(x.network, indicator_node, param_key)[mat]
    #loop through orders on relevant arc
    if src == dst || early_fulfillment #assume production at plant accepts early fulfillment
        rel_orders = all_relevant_orders(x,src,dst,mat)
    else #only attempt to fulfill orders that have expired their service lead time
        rel_orders = expired_relevant_orders(x,src,dst,mat)
    end
    #return sorted by soonest due date
    return sort(rel_orders, [:due, :id], view=true) #sort by service lead time and creation date (serve orders with lowest service lead time, ranked by order age)
end
"""
    all_relevant_orders(x::SupplyChainEnv, src::Int, [dst::Union{Int,Symbol},] mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` with any due date.
"""
all_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}) =
    @rsubset(x.open_orders, :material == mat, :arc == (src,dst), view=true)

"""
    expired_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` that are due now or expired.
"""
expired_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}) = 
    @rsubset(x.open_orders, :due <= 0, :material == mat, :arc == (src,dst), view=true)

"""
    raw_material_orders(x::SupplyChainEnv, ids, rmat_names)

Get a grouped dataframe of raw material consumption orders (:consumption) associated with order ids in `ids`
"""
raw_material_orders(x::SupplyChainEnv, ids, rmat_names) =
    @chain x.open_orders begin
        @rsubset(:id in Set(ids), :material in Set(rmat_names), view=true)
        groupby([:id,:material])
    end

"""
    reallocate_demand!(x::SupplyChainEnv, order_row::DataFrameRow)

Calculate unfulfilled demand and reallocate order to next supplier in the priority list.
"""
function reallocate_demand!(x::SupplyChainEnv, order_row::DataFrameRow)
    sup, req = order_row.arc
    if order_row.amount <= 0 || order_row.due > 0 || sup == req || req == :market #if order amount is not positive OR order is not due OR production order OR market demand, do NOT reallocate
        return
    end
    mat = order_row.material
    new_alloc = missing
    sup_priority = get_prop(x.network, req, :supplier_priority)[mat]
    if length(sup_priority) > 1
        sup_priority = vcat(sup_priority, sup_priority[1]) #add top priority to the end of the list so that lowest priority reallocates to top priority (cycle)
        next_sup_loc = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
        new_sup = sup_priority[next_sup_loc] #next supplier in line
        new_alloc = (new_sup, req) #store new arc where reallocated
        order_row.arc = new_alloc #update x.open_orders to reallocate material
    end
end

"""
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    orders_grp = groupby(copy(x.open_orders), [:arc, :material])

    for n in x.markets
        dmnd_freq_dict = get_prop(x.network, n, :demand_frequency)
        dmnd_dist_dict = get_prop(x.network, n, :demand_distribution)
        serv_dist_dict = get_prop(x.network, n, :service_lead_time)
        n_mats = get_prop(x.network, n, :node_materials)
        for mat in n_mats
            #get demand parameters
            p = dmnd_freq_dict[mat] #probability of demand occuring
            dmnd = dmnd_dist_dict[mat] #demand distribution
            serv_lt = serv_dist_dict[mat] #service lead time
            last_order_id = x.num_orders #last order created in system
            #place p orders
            for _ in 1:floor(p) + rand(Bernoulli(p % 1)) #p is the number of orders. If fractional, the fraction is the likelihood of rounding up.
                q = rand(dmnd) #sample demand
                slt = ceil(Int,rand(serv_lt)) #sample service lead time
                if q > 0
                    create_order!(x, n, :market, mat, q, slt)
                end
            end
            #if no new orders were created, check if there are any pending orders; if not, the exit iteration
            if last_order_id == x.num_orders
                pending = calculate_pending([(n,:market)], mat, orders_grp)
                iszero(pending) && continue
            end
            #try to fulfill demand from stock (will include any open orders)
            fulfill_from_stock!(x, n, :market, mat, 0) #0 lead time since at market
            #if MTO and lead time is 0, try to fulfill from production
            if ismto(x,n,mat) #fulfill make-to-order from production (NOTE: COULD RECHECK IF THERE IS ANY PENDING ORDER THAT WASN"T FULFILLED FROM STOCK)
                lt_dist = get_prop(x.network, n, n, :lead_time)[mat]
                if iszero(lt_dist)
                    capacities = get_prop(x.network, n, :production_capacity)
                    prod_lt = ceil(Int,rand(lt_dist))
                    fulfill_from_production!(x, n, :market, mat, prod_lt, capacities)
                end
            end
        end
    end

    #lost sales for expired orders (if backlog = false, or guaranteed service = true)
    if !x.options[:backlog] || x.options[:guaranteed_service] #any open orders with non-positive relative due date (0 service lead time) become lost sales
        lost_sales!(x)
    end
end