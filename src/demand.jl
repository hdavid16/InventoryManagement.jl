"""
    place_orders!(x::SupplyChainEnv, act::NamedArray)

Place inventory replenishment requests throughout the network.
"""
function place_orders!(x::SupplyChainEnv, act::NamedArray)
    arcs = names(act,2) #arcs
    #exit if no action and there is no backlog
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

    #place requests
    for req in request_nodes #loop by nodes placing requests (in reverse topological order)
        sup_priority = get_prop(x.network, req, :supplier_priority) #get supplier priority list
        for mat in mats, sup in sup_priority[mat] #loop by material and supplier priority
            a = (sup, req) #arc
            lead = float(leads[Edge(a...), mat]) #sampled lead time
            serv = float(servs[Edge(a...), mat]) #sampled service lead time
            amount = act[mat, a] #amount requested
            #continue to next iteration if no request made and no backlog
            if iszero(amount)
                #calculate outstanding orders (backlog) to determine if an order should be placed even if amount = 0 (will be 0 if assuming lost sales)
                backlog = calculate_backlog(x, a..., mat)
                if req in x.markets
                    backlog += calculate_backlog(x, req, [:market], mat)
                end
                if req in x.producers
                    backlog += calculate_backlog(x, req, [:production], mat)
                end
                if iszero(backlog)
                    push!(x.demand, [x.period, a, mat, 0, 0, 0, 0, missing])
                    continue 
                end
            end
            #create order and save service lead time
            create_order!(x, a..., mat, amount, serv)
            if sup == req
                fulfill_from_production!(x, a..., mat, lead, supply_grp, pipeline_grp, capacities)
            else
                fulfill_from_stock!(x, a..., mat, lead, supply_grp, pipeline_grp)
            end
            #check for any backlogged market demand orders
            if x.options[:backlog] && req in x.markets 
                fulfill_from_stock!(x, req, :market, mat, 0., supply_grp, missing)
            end
        end
    end

    #lost sales for expired orders (if backlog = false, or guaranteed service = true)
    if !x.options[:backlog] || x.options[:guaranteed_service] #any open orders with non-positive relative due date (0 service lead time) become lost sales
        expired_orders = filter(:due => t -> t <= 0, x.open_orders, view=true).id
        transform!(
            filter(:id => id -> id in expired_orders, x.orders, view=true),
            [:arc, :fulfilled] => ByRow((a,log) -> vcat(log, (time=x.period, supplier=a[1], amount=:lost_sale))) => :fulfilled
        )
        filter!(:due => t -> t > 0, x.open_orders)
    end

    #updated production capacities (after commited production)
    for n in x.producers
        set_prop!(x.network, n, :production_capacity, capacities[n])
    end
end

"""
    exit_place_orders!(x, arcs)

Abort order placement.
"""
function exit_place_orders!(x::SupplyChainEnv, arcs::Vector)
    for a in arcs, mat in x.materials
        backlog = calculate_backlog(x, a..., mat)
        push!(x.demand, [x.period, a, mat, 0, 0, 0, backlog, missing])
    end
end

"""
    create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Float64)

Create and log order.
"""
function create_order!(x::SupplyChainEnv, sup::Int, req::Union{Int,Symbol}, mat::Union{Symbol,String}, amount::Real, service_lead_time::Float64)
    x.num_orders += 1 #create new order ID
    push!(x.orders, [x.num_orders, x.period, (sup,req), mat, amount, []]) #update order history
    push!(x.open_orders, [x.num_orders, (sup,req), mat, amount, service_lead_time]) #add order to temp order df
    if sup == req && isproduced(x, sup, mat)
        bom = get_prop(x.network, sup, :bill_of_materials)
        rmat_names = names(filter(k -> k < 0, bom[:,mat]), 1) #names of raw materials
        for rmat in rmat_names
            push!(x.open_orders, [x.num_orders, (sup,:production), rmat, -amount*bom[rmat,mat], service_lead_time])
        end
    end
    sort!(x.open_orders, [:due, :id]) #sort by service lead time and creation date (serve orders with lowest service lead time, ranked by order age)
end

"""
    relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})

Return filtered dataframe of active orders relevant to the arc `(src, dst)` for material `mat`
"""
function relevant_orders(x::SupplyChainEnv, src::Int, dst::Union{Int,Symbol}, mat::Union{Symbol,String})
    #node from which to check early_fulfillment param value
    indicator_node = dst == :market ? src : dst
    #loop through orders on relevant arc
    if get_prop(x.network, indicator_node, :early_fulfillment)[mat]
        orders_df = filter(
            [:arc, :material] => 
                (a,m) -> 
                    a == (src,dst) && #nodes with the arc
                    m == mat, #orders for the same material
            x.open_orders, 
            view = true
        )
    else #only attempt to fulfill orders that have expired their service lead time
        orders_df = filter(
            [:arc, :material, :due] => 
                (a,m,t) -> 
                    a == (src,dst) && #nodes with the arc
                    m == mat &&  #orders for the same material
                    t <= 0, #expired orders
            x.open_orders, 
            view = true
        )
    end

    return orders_df
end

"""
    log_unfulfilled_demand!(x::SupplyChainEnv, order_row::DataFrameRow, accepted::Float64)

Calculate unfulfilled demand and reallocate order to next supplier in the priority list.
"""
function log_unfulfilled_demand!(x::SupplyChainEnv, order_row::DataFrameRow, accepted::Float64)
    sup, req = order_row.arc
    mat = order_row.material
    new_alloc = missing
    if x.options[:reallocate] && (sup != req || req != :market) #don't reallocate if external demand or production request
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
    for n in x.markets
        dmnd_freq_dict = get_prop(x.network, n, :demand_frequency)
        dmnd_dist_dict = get_prop(x.network, n, :demand_distribution)
        serv_dist_dict = get_prop(x.network, n, :service_lead_time)
        for mat in x.materials
            p = dmnd_freq_dict[mat] #probability of demand occuring
            dmnd = dmnd_dist_dict[mat] #demand distribution
            serv_lt = serv_dist_dict[mat] #service lead time
            #place p orders
            for _ in 1:floor(p) + rand(Bernoulli(p % 1)) #p is the number of orders. If fractional, the fraction is the likelihood of rounding up.
                q = rand(dmnd) #sample demand
                slt = rand(serv_lt) #sample service lead time
                if q > 0
                    external_order!(x,n,mat,q,slt)
                end
            end
        end
    end
end

"""
    external_order!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, q::Real, serv::Real)

Create external demand at node `n` for material `mat` for quantity `q` with service lead time `serv`.
"""
function external_order!(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, q::Real, serv::Real)
    supply_df = filter([:period,:node] => (t,n) -> t == x.period && n in x.markets, x.inventory_on_hand, view=true) #on_hand inventory at node
    supply_grp = groupby(supply_df, [:node, :material]) #group by node and material
    create_order!(x, n, :market, mat, q, serv)
    fulfill_from_stock!(x, n, :market, mat, 0., supply_grp, missing) #0 lead time since at market; pipeline_grp is missing since no arc betwen market node and market
end