"""
    replenishment_orders!(x::SupplyChainEnv, act::NamedArray)

Place inventory replenishment requests throughout the network.
"""
function replenishment_orders!(x::SupplyChainEnv, act::NamedArray)
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
    supplier_nodes = filter( #nodes supplying inventory
        n -> !isempty(outneighbors(x.network, n)), #exclude sink nodes
        topological_sort(x.network) #sort nodes in topological order so that orders are placed moving down the network
    )
    #get copy of existing open orders
    orders_grp = groupby(copy(x.open_orders), [:arc, :material])

    #place requests
    while !isempty(supplier_nodes)
        sup = popfirst!(supplier_nodes)
        reqs = outneighbors(x.network, sup) #requestors (downstream nodes). 
        arcs = [(sup,req) for req in reqs]
        for mat in get_prop(x.network, sup, :node_materials) #loop by material
            #continue to next iteration if no request made and no backlog
            if iszero(sum(act[mat,arcs]))
                #determine if there are any open orders for this material on the arc. If not, skip the fulfillment step.
                pending = calculate_pending(
                    vcat(arcs,(sup,:consumption)), #include any raw material conversion orders (:consumption)
                    mat, 
                    orders_grp
                )
                if iszero(pending)
                    continue 
                end
            else
                #create order and save service lead time
                for arc in arcs
                    amount = act[mat, arc] #amount requested
                    if amount > 0
                        serv = servs[Edge(arc...), mat] #sampled service lead time
                        create_order!(x, arc..., mat, amount, serv)
                    end
                end
            end
            #fulfill production orders
            if sup in x.producers && isproduced(x.network, sup, mat)
                lead = leads[Edge(sup,sup), mat] #sampled lead time
                fulfill_from_production!(x, sup, sup, mat, lead, capacities[sup])
            end 
            #fulfill stock transfer orders
            reallocated_suppliers = fulfill_from_stock!(x, sup, setdiff(reqs,sup), mat, leads)
            pushfirst!(supplier_nodes, reallocated_suppliers...) #preappend reallocated suppliers
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
    expired_orders = filter(:due => <=(0), x.open_orders, view=true)
    for row in eachrow(expired_orders)
        push!(x.fulfillments, (row.id, x.period, row.arc, row.material, row.amount, :lost_sale))
    end
    if !isempty(expired_orders) #remove expired orders if there are any
        lost_sales = Set(expired_orders.id)
        filter(:id => in(lost_sales), x.orders, view=true).fulfilled .= :lost_sale #log lost sale
        filter!(:due => >(0), x.open_orders)
    end
end

"""
    create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Material, amount::Real, service_lead_time::Float64)

Create and log order.
"""
function create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Material, amount::Real, service_lead_time::Real)
    x.num_orders += 1 #create new order ID
    push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,req), mat, amount, missing]) #update order history
    push!(x.open_orders, [x.num_orders, x.period, service_lead_time, (sup,req), mat, amount]) #add order to temp order df
    #create raw material consumption and coproduction orders
    if (sup == req && isproduced(x,sup,mat)) || (req == :market && ismto(x,sup,mat))
        bom = get_prop(x.network, sup, :bill_of_materials)
        rmat_names = names(filter(<(0), bom[:,mat]), 1) #names of raw materials
        cmat_names = names(filter(>(0), bom[:,mat]), 1) #names of raw materials
        for rmat in rmat_names
            push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,:consumption), rmat, -amount*bom[rmat,mat], missing]) #update order history
            push!(x.open_orders, [x.num_orders, x.period, service_lead_time, (sup,:consumption), rmat, -amount*bom[rmat,mat]])
        end
        for cmat in cmat_names
            push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,:coproduction), cmat, amount*bom[cmat,mat], missing]) #update order history
            push!(x.open_orders, [x.num_orders, x.period, service_lead_time, (sup,:coproduction), cmat, amount*bom[cmat,mat]])
        end
    end
end

"""
    relevant_orders(x::SupplyChainEnv, src::Int, dsts::Vector{Int}, mat::Material)

Return filtered dataframe of active orders relevant to the arcs downstream of `src` for material `mat`.

    relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material)

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat`.
"""
function relevant_orders(x::SupplyChainEnv, src::Int, dsts::Vector{Int}, mat::Material)
    arcs = Set([(src,dst) for dst in dsts]) #build all arcs
    #filter by arc and due date if early fulfillment not allowed (if early fulfillment not allowed, only include due/expired orders)
    mask(arc,due) = arc in arcs && (get_prop(x.network, arc[2], :early_fulfillment)[mat] || due <= 0)
    rel_orders = subset(x.open_orders,
        :material => ByRow(==(mat)),
        :arc => ByRow(in(arcs)),
        [:arc, :due] => ByRow(mask),
        view=true
    )
    customer_priority = get_prop(x.network, src, :customer_priority)[mat]
    sort( #return sorted by soonest due date
        rel_orders, #NOTE: Can add customer priority when sorting
        [
            order(:due, by = t -> iszero(t) ? -Inf : t), #give highest priority to order due now (due = 0), then rank by due date (expired orders if any are ranked higher), 
            :created, #then rank by order creation date
            order(:arc, by = a -> findfirst(==(a[2]), customer_priority)), #prioritize by customer priority
        ], #NOTE: ASSUME SORT FIRST BY TIME AND THEN BY PRIORITY!
        view=true
    ) 
end
function relevant_orders(x::SupplyChainEnv, src::Int, dst::Int, mat::Material)
    rel_orders = all_relevant_orders(x,src,dst,mat) #assume production at plant accepts early fulfillment
    sort( #return sorted by soonest due date
        rel_orders, 
        [order(:due, by = t -> iszero(t) ? -Inf : t), :created], #give highest priority to order due now (due = 0), then rank by due date (expired orders if any are ranked higher), then rank by order creation date
        view=true
    ) 
end
function relevant_orders(x::SupplyChainEnv, src::Int, dst::Symbol, mat::Material)
    if get_prop(x.network, src, :market_early_fulfillment)[mat]
        rel_orders = all_relevant_orders(x,src,dst,mat)
    else #only attempt to fulfill orders that have expired their service lead time
        rel_orders = due_relevant_orders(x,src,dst,mat)
    end
    sort( #return sorted by soonest due date
        rel_orders, 
        [order(:due, by = t -> iszero(t) ? -Inf : t), :created], #give highest priority to order due now (due = 0), then rank by due date (expired orders if any are ranked higher), then rank by order creation date
        view=true
    ) 
end
"""
    all_relevant_orders(x::SupplyChainEnv, src::Int, [dst::Union{Int,Symbol},] mat::Material)

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` with any due date.
"""
function all_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material)
    a = (src,dst)
    subset(x.open_orders,
        :material => ByRow(==(mat)),
        :arc => ByRow(==(a)),
        view=true
    )
end

"""
    due_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material)

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` that are due now or expired.
"""
function due_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Material)
    a = (src,dst)
    subset(x.open_orders,
        :due => ByRow(<=(0)),
        :material => ByRow(==(mat)),
        :arc => ByRow(==(a)),
        view=true
    )
end

"""
    raw_material_orders(x::SupplyChainEnv, ids, rmat_names)

Get a grouped dataframe of raw material consumption orders (:consumption) associated with order ids in `ids`
"""
function raw_material_orders(x::SupplyChainEnv, ids, rmat_names)
    id_set = Set(ids)
    rmat_set = Set(rmat_names)
    groupby(
        subset(x.open_orders,
            :id => ByRow(in(id_set)),
            :material => ByRow(in(rmat_set)),
            view=true
        ),
        [:id, :material]
    )
end

"""
    update_orders!(x::SupplyChainEnv, orders_df::SubDataFrame)

Update system orders (remove fulfilled orders from `x.open_orders` and log order fulfillment in `x.orders`).
"""
function update_orders!(x::SupplyChainEnv, orders_df::SubDataFrame)
    fulfilled_order_list = filter(:amount => <=(0), orders_df, view=true).id
    if !isempty(fulfilled_order_list)
        fulfilled_orders = Set(fulfilled_order_list) #find fulfilled orders
        if length(fulfilled_order_list) != length(fulfilled_orders) #there are production orders that were fulfilled from stock (delete these)
            prod_keys = Set([:consumption,:coproduction])
            filter!([:id,:arc] => (id,arc) -> !(id in fulfilled_orders && arc[2] in prod_keys), x.orders)
        end
        filter(:id => in(fulfilled_orders), x.orders, view=true).fulfilled .= x.period #store fulfillment date
        filter!(:id => !in(fulfilled_orders), x.open_orders) #remove fulfilled orders
    end
end

"""
    reallocate_demand!(x::SupplyChainEnv, order_row::DataFrameRow)

Calculate unfulfilled demand and reallocate order to next supplier in the priority list.
"""
function reallocate_demand!(x::SupplyChainEnv, order_row::DataFrameRow)
    if order_row.amount <= 0 || order_row.due > 0 #if order amount is not positive OR order is not due, do NOT reallocate
        return
    end
    sup, req = order_row.arc
    mat = order_row.material
    new_sup = nothing
    sup_priority = get_prop(x.network, req, :supplier_priority)[mat]
    if length(sup_priority) > 1
        next_sup_loc = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
        if next_sup_loc <= length(sup_priority) #stop reallocating when you have reached the lowest priority supplier
            new_sup = sup_priority[next_sup_loc] #next supplier in line
            new_alloc = (new_sup, req) #store new arc where reallocated
            order_row.arc = new_alloc #update x.open_orders to reallocate material
        end
    end

    return new_sup
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
            if x.tmp[n,mat,:on_hand] > 0
                fulfill_from_stock!(x, n, :market, mat) #0 lead time since at market
            end
            #if MTO and lead time is 0, try to fulfill from production
            if last_order_id > x.num_orders && ismto(x,n,mat) #fulfill make-to-order from production if at least 1 new MTO created (NOTE: COULD RECHECK IF THERE IS ANY PENDING ORDER THAT WASN"T FULFILLED FROM STOCK)
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