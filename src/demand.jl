"""
    place_orders!(x::SupplyChainEnv, act::NamedArray)

Place inventory replenishment requests throughout the network.
"""
function place_orders!(x::SupplyChainEnv, act::NamedArray)
    arcs = names(act,2) #arcs

    #exit if no action and there is no backlogging
    if !x.options[:backlog] && iszero(act)
        exit_place_orders!(x, arcs)
        return
    end
    #extract info
    mats = x.materials #materials
    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    #sample lead times and service times
    leads = Dict((a,mat) => rand(get_prop(x.network, a, :lead_time)[mat]) for a in edges(x.network), mat in mats)
    servs = Dict((a,mat) => rand(get_prop(x.network, a, :service_lead_time)[mat]) for a in edges(x.network), mat in mats)
    #identify nodes that can place requests
    nodes = topological_sort(x.network) #sort nodes in topological order so that orders are placed moving down the network
    source_nodes = filter(n -> isempty(inneighbors(x.network, n)), nodes) #source nodes (can't place replenishment orders)
    request_nodes = setdiff(nodes, source_nodes) #nodes placing requests (all non-source nodes)
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inventory_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inventory_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])
    #get existing open orders
    orders_grp = groupby(copy(x.open_orders), [:arc, :material])

    #place requests
    for req in request_nodes #loop by nodes placing requests (in reverse topological order)
        sup_priority = get_prop(x.network, req, :supplier_priority) #get supplier priority list
        req_mats = get_prop(x.network, req, :node_materials) #get materials that are compatible with requester node
        for mat in req_mats, sup in sup_priority[mat] #loop by material and supplier priority
            a = (sup, req) #arc
            lead = float(leads[Edge(a...), mat]) #sampled lead time
            serv = float(servs[Edge(a...), mat]) #sampled service lead time
            amount = act[mat, a] #amount requested
            #continue to next iteration if no request made and no backlog
            if iszero(amount)
                #calculate outstanding orders (backlog) to determine if an order should be placed even if amount = 0 (will be 0 if assuming lost sales)
                arcs_list = [a,(req,:production)]
                backlog = calculate_backlog(x, arcs_list, mat, orders_grp) #include any outstanding raw material conversion orders (:production)
                if iszero(backlog)
                    push!(x.demand, [x.period, a, mat, 0, 0, 0, 0, missing])
                    continue 
                end
            else
                #create order and save service lead time
                create_order!(x, a..., mat, amount, serv)
            end
            #try to fulfill request
            if sup == req
                fulfill_from_production!(x, a..., mat, lead, supply_grp, pipeline_grp, capacities[sup])
            else
                fulfill_from_stock!(x, a..., mat, lead, supply_grp, pipeline_grp)
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
    expired_orders = filter(:due => t -> t <= 0, x.open_orders, view=true)
    for row in eachrow(expired_orders)
        push!(x.fulfillments, (row.id, x.period, row.arc[1], :lost_sale))
    end
    if !isempty(expired_orders)
        filter!(:due => t -> t > 0, x.open_orders)
    end
end

"""
    exit_place_orders!(x::SupplyChainEnv, arcs::Vector)

Abort order placement.
"""
function exit_place_orders!(x::SupplyChainEnv, arcs::Vector)
    due_by = x.options[:adjusted_stock] ? Inf : 0 #Inf means that all orders placed are counted (even if not due); otherwise, only due orders are counted
    orders_df = filter(:due => j -> j <= due_by, x.open_orders, view=true) #orders to count in backlogging
    orders_grp = groupby(orders_df, [:arc, :material]) #group by arc and material
    for a in arcs
        dst_mats = get_prop(x.network, a[2], :node_materials)
        for mat in dst_mats
            backlog = calculate_backlog(x, [a], mat, orders_grp)
            push!(x.demand, [x.period, a, mat, 0, 0, 0, backlog, missing])
        end
    end
end

"""
    create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Float64)

Create and log order.
"""
function create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Real)
    x.num_orders += 1 #create new order ID
    push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,req), mat, amount]) #update order history
    push!(x.open_orders, [x.num_orders, (sup,req), mat, amount, service_lead_time]) #add order to temp order df
    if (sup == req && isproduced(x,sup,mat)) || (req == :market && ismto(x,sup,mat))
        bom = get_prop(x.network, sup, :bill_of_materials)
        rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
        for rmat in rmat_names
            push!(x.orders, [x.num_orders, x.period, x.period + service_lead_time, (sup,:production), rmat, -amount*bom[rmat,mat]]) #update order history
            push!(x.open_orders, [x.num_orders, (sup,:production), rmat, -amount*bom[rmat,mat], service_lead_time])
        end
    end
    sort!(x.open_orders, [:due, :id]) #sort by service lead time and creation date (serve orders with lowest service lead time, ranked by order age)
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
        return all_relevant_orders(x,src,dst,mat)
    else #only attempt to fulfill orders that have expired their service lead time
        return expired_relevant_orders(x,src,dst,mat)
    end
end
"""
    all_relevant_orders(x::SupplyChainEnv, src::Int, [dst::Union{Int,Symbol},] mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` with any due date.
"""
all_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}) =
    filter( 
        [:arc, :material] => 
            (a,m) -> 
                a == (src,dst) && #production orders at `n`
                m == mat, #orders for the same material
        x.open_orders, 
        view = true
    )
"""
    expired_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat` that are due now or expired.
"""
expired_relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String}) = 
    filter(
        [:arc, :material, :due] => 
            (a,m,t) -> 
                a == (src,dst) && #nodes with the arc
                m == mat &&  #orders for the same material
                t <= 0, #expired orders
        x.open_orders, 
        view = true
    )

"""
    raw_material_orders(x::SupplyChainEnv, ids, rmat_names)

Get a grouped dataframe of raw material consumption orders (:production) associated with order ids in `ids`
"""
function raw_material_orders(x::SupplyChainEnv, ids, rmat_names)
    raw_orders_df = filter(
        [:id, :material] => 
            (id,m) -> 
                id in ids && 
                m in rmat_names, 
        x.open_orders, 
        view = true
    )
    return groupby(raw_orders_df, [:id, :material])
end

"""
    log_unfulfilled_demand!(x::SupplyChainEnv, order_row::DataFrameRow, accepted::Float64)

Calculate unfulfilled demand and reallocate order to next supplier in the priority list.
"""
function log_unfulfilled_demand!(x::SupplyChainEnv, order_row::DataFrameRow, accepted::Float64)
    sup, req = order_row.arc
    mat = order_row.material
    new_alloc = missing
    if x.options[:reallocate] && sup != req && req != :market #don't reallocate if external demand or production request
        sup_priority = get_prop(x.network, req, :supplier_priority)[mat]
        if length(sup_priority) > 1
            sup_priority = vcat(sup_priority, sup_priority[1]) #add top priority to the end of the list so that lowest priority reallocates to top priority (cycle)
            next_sup_loc = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
            new_sup = sup_priority[next_sup_loc] #next supplier in line
            new_alloc = (new_sup, req) #store new arc where reallocated
            order_row.arc = new_alloc #update x.open_orders to reallocate material
        end
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
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    supply_df = filter([:period,:node] => (t,n) -> t == x.period && n in x.markets, x.inventory_on_hand, view=true) #on_hand inventory at node
    supply_grp = groupby(supply_df, [:node, :material]) #group by node and material

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
            #place p orders
            for _ in 1:floor(p) + rand(Bernoulli(p % 1)) #p is the number of orders. If fractional, the fraction is the likelihood of rounding up.
                q = rand(dmnd) #sample demand
                slt = rand(serv_lt) #sample service lead time
                if q > 0
                    create_order!(x, n, :market, mat, q, slt)
                end
            end
            #fulfill demand (will include any backlogged orders)
            #if material is make-to-order and the production time is 0, then trigger production directly
            trigger_production = 
                ismto(x,n,mat) && #check if material is make-to-order
                (get_prop(x.network, n, n, :lead_time)[mat] |> 
                    lt -> !isa(lt, Sampleable) && iszero(lt)) #0 production lead time
            if trigger_production
                capacities = get_prop(x.network, n, :production_capacity)
                fulfill_from_production!(x, n, :market, mat, 0., supply_grp, missing, capacities)
            #fulfill orders from stock 
            else
                fulfill_from_stock!(x, n, :market, mat, 0., supply_grp, missing) #0 lead time since at market; pipeline_grp is missing since no arc betwen market node and market
            end
        end
    end

    #lost sales for expired orders (if backlog = false, or guaranteed service = true)
    if !x.options[:backlog] || x.options[:guaranteed_service] #any open orders with non-positive relative due date (0 service lead time) become lost sales
        lost_sales!(x)
    end
end