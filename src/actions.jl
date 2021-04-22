function show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)
    prods = env.products
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    act = reshape(action, (length(prods), length(arcs)))
    df = DataFrame(:product => prods,
                   [Symbol(a) => act[:,i] for (i,a) in enumerate(arcs)]...)
end

function (x::SupplyChainEnv)(action::Vector{T} where T <: Real)
    #validate action input
    @assert all(action .>= 0) "Reorder actions cannot be negative."
    @assert length(action) == length(x.products)*ne(x.network) "Reorder action vector must have length num_products * num_edges."

    #increase period counter
    x.period += 1

    #intialize next period on-hand and pipeline inventories with previous inventories
    for n in union(x.markets, x.distributors), p in x.products
        prev = filter(j -> j.period == x.period-1 && j.node == n && j.product == p, x.inv_on_hand).level[1] #previous inventory level
        push!(x.inv_on_hand, [x.period, n, p, prev]) #intialize with previous inventory levels
    end
    for a in edges(x.network), p in x.products
        prev = filter(j -> j.period == x.period-1 && j.arc == (a.src,a.dst) && j.product == p, x.inv_pipeline).level[1] #previous inventory level
        push!(x.inv_pipeline, [x.period, (a.src,a.dst), p, prev]) #intialize with previous inventory levels
    end

    #move active shipments forward one period
    x.shipments.lead .-= 1

    #move active production forward one period
    x.production.lead .-= 1

    #sample lead times
    leads = Dict(a => rand(get_prop(x.network, a, :lead_time)) for a in edges(x.network)) #sample lead time distribution

    #reshape action
    prods = x.products
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    act = reshape(action, (length(prods), length(arcs)))

    #store original production capacities and inventory levels (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    inv_levels = Dict(n => filter(i -> i.period == x.period && i.node == n, x.inv_on_hand)[:,[:product, :level]]
                        for n in x.distributors)

    #place requests
    for i in 1:length(prods) #loop by products
        p = prods[i] #product requested
        for n in union(x.distributors, x.markets) #loop by nodes placing requests
            sup_priority = get_prop(x.network, n, :supplier_priority)[p] #get supplier priority list
            for sup in sup_priority #loop by supplier priority
                a = (sup, n) #arc
                j = findfirst(k -> k == a, arcs) #find index for that arc in the action matrix
                amount = act[i,j] #amount requested
                prod_time = 0 #initialize production time

                #accept or adjust requests
                if a[1] in x.producers
                    capacity = capacities[a[1]][p] #get production capacity
                    accepted = min(amount, capacity) #accepted request
                    act[i,j] = accepted #update order quantity
                    capacities[a[1]][p] = capacity - accepted #update production capacity to account for commited capacity (handled first come first serve)
                    if accepted > 0 #schedule production
                        prod_time = get_prop(x.network, a[1], :production_time)[p] #production time
                        push!(x.production, (a, p, accepted, prod_time))
                    end
                else
                    supply = filter(i -> i.product == p, inv_levels[a[1]]).level[1] #on_hand inventory at supplier
                    accepted = min(amount, supply) #accepted request
                    act[i,j] = accepted #update order quantity
                    inv_levels[a[1]][inv_levels[a[1]].product .== p, :level] .= supply - accepted #update available inventory to account for commited inventory (handled first come first serve)
                    if accepted > 0 #subtract sent onhand inventory, add sent pipeline inventory
                        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                      (x.inv_on_hand.node .== a[1]) .&
                                      (x.inv_on_hand.product .== p), :level] .-= accepted
                        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                                       (x.inv_pipeline.product .== p), :level] .+= accepted
                    end
                end

                #chekc if some was not accepted and reallocate
                unfulfilled = amount - accepted #amount not accepted
                if unfulfilled > 0 #reallocate unfulfilled request to next priority supplier
                    @warn "Replenishment request made by node $(a[2]) to node $(a[1]) for product $p was reduced from $amount to $accepted due to insufficient production capacity or insufficient inventory."
                    if x.reallocate
                        next_sup = findfirst(k -> k == a[1], sup_priority) + 1 #get next in line
                        if next_sup <= length(sup_priority) #check that there is a next one in line
                            @warn "Reallocating non-accepted request for node $(a[2]) for product $p to node $next_sup (amount = $unfulfilled)"
                            jj = findfirst(k -> k == (next_sup, a[2]), arcs) #find index for that arc in the action matrix
                            act[i,jj] += unfulfilled #add unfulfilled to the next supplier in the line
                        end
                    end
                end

                #store shipments and lead times
                if accepted > 0
                    lead = leads[Edge(a[1],a[end])] + prod_time
                    push!(x.shipments, [a, p, accepted, lead]) #update active shipments
                else
                    lead = 0 #no request made
                end
                push!(x.replenishments, [x.period, a, p, accepted, lead])
            end
        end
    end

    #update production and associated pipeline inventories
    produced = filter(i -> iszero(i.lead), x.production) #find active production that has completed
    for i in 1:nrow(produced) #add produced material to pipeline
        a, p, amount = produced[i, [:arc, :product, :amount]]
        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                       (x.inv_pipeline.product .== p), :level] .+= amount
    end
    filter!(i -> i.lead > 0, x.production) #remove produced inventory that was shipped

    #update on hand and pipeline inventories due to arrived shipments
    arrivals = filter(i -> iszero(i.lead), x.shipments) #find active shipments with 0 lead time
    for i in 1:nrow(arrivals)
        a, p, amount = arrivals[i, [:arc, :product, :amount]]
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== a[end]) .&
                      (x.inv_on_hand.product .== p), :level] .+= amount
        x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                       (string.(x.inv_pipeline.arc) .== string(a)) .&
                       (x.inv_pipeline.product .== p), :level] .-= amount
    end
    filter!(i -> i.lead > 0, x.shipments) #remove shipments that arrived

    #update inventory positions for non-market nodes
    for n in x.distributors, p in x.products #update distribution centers (not plants, not end distributors [markets])
        upstream = sum(filter(j -> j.period == x.period && j.arc[end] == n && j.product == p, x.inv_pipeline).level)
        onhand = filter(j -> j.period == x.period && j.node == n && j.product == p, x.inv_on_hand).level[1]
        push!(x.inv_position, [x.period, n, p, onhand + upstream]) #update inventory
    end

    #markets open and demand occurs
    for m in x.markets, p in x.products
        freq = Bernoulli(get_prop(x.network, m, :demand_frequency)[p]) #demand frequency
        dmnd = get_prop(x.network, m, :demand_distribution)[p] #demand distribution
        q = rand(freq) * rand(dmnd) #quantity requested
        demand = [x.period, m, p, q, 0, 0, 0] #initialize demand vector to store in df
        inv = filter(i -> i.period == x.period && i.node == m && i.product == p, x.inv_on_hand).level[1] #current inventory at node
        if x.backlog #add previous backlog to the quantity requested at the market
            q += filter(i -> i.period == x.period-1 && i.node == m && i.product == p, x.demand).backlog[1]
        end
        demand[5] = min(q, inv) #sales
        demand[6] = max(q - inv, 0) #unfilfilled
        demand[7] = x.backlog ? demand[6] : 0 #backlog
        push!(x.demand, demand) #update df

        #update end of period inventory (subtract sales)
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== m) .&
                      (x.inv_on_hand.product .== p), :level] .-= demand[5]

        #update inventory position
        upstream = sum(filter(j -> j.period == x.period && j.arc[end] == m && j.product == p, x.inv_pipeline).level)
        onhand = filter(j -> j.period == x.period && j.node == m && j.product == p, x.inv_on_hand).level[1]
        push!(x.inv_position, [x.period, m, p, onhand + upstream - demand[7]]) #update inventory (include backlog)
    end

    #calculate profit at each node
    for n in vertices(x.network)
        profit = 0 #initialize node profit
        for p in x.products
            #get costs
            if !(n in x.producers) #holding cost
                hold_cost = get_prop(x.network, n, :holding_cost)[p]
                onhand = filter(j -> j.period == x.period && j.node == n && j.product == p, x.inv_on_hand).level[1]
                profit -= hold_cost * onhand
            end
            if n in x.producers #production cost
                prod_cost = get_prop(x.network, n, :production_cost)[p]
                produced = sum(filter(j -> j.period == x.period && j.arc[1] == n && j.product == p, x.replenishments).amount)
                profit -= prod_cost * produced
            elseif n in x.markets #sales profit at end distributors (and penalize for unfulfilled demand)
                sales_price = get_prop(x.network, n, :sales_price)[p]
                dmnd_penalty = get_prop(x.network, n, :demand_penalty)[p]
                sold, unfilled = filter(j -> j.period == x.period && j.node == n && j.product == p, x.demand)[1, [:sale, :unfulfilled]]
                profit += sales_price * sold
                profit -= dmnd_penalty * unfilled
            end
            #pay suppliers for received inventory and pay transportation cost
            for pred in inneighbors(x.network, n)
                price = get_prop(x.network, pred, n, :sales_price)[p]
                trans_cost = get_prop(x.network, pred, n, :transportation_cost)[p]
                #pay purchase of inventory
                purchased = filter(j -> j.arc == (pred, n) && j.product == p, arrivals).amount
                if !isempty(purchased)
                    profit -= purchased[1] * price
                end
                #pay transportation (pipeline hoding) cost
                intransit = filter(j -> j.period == x.period && j.arc == (pred, n) && j.product == p, x.inv_pipeline).level[1]
                profit -= intransit * trans_cost
            end
            #receive payment for delivered inventory and pay production costs
            for succ in outneighbors(x.network, n)
                price = get_prop(x.network, n, succ, :sales_price)[p]
                #receive payment for delivered inventory
                sold = filter(j -> j.arc == (n, succ) && j.product == p, arrivals).amount
                if !isempty(sold)
                    profit += sold[1] * price
                end
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end

    #update reward (current profit). Should be revised for RL
    x.reward = sum(filter(j -> j.period == x.period, x.profit).value)
end

select_any_supplier(::T) where T = true

function reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                level::Symbol = :position, kind::Symbol = :rQ,
                supplier_selection::Symbol = :priority)

    #read parameters
    t = env.period
    nodes = union(env.distributors, env.markets)
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    prods = env.products

    #check inputs
    @assert kind in [:rQ, :sS] "The policy kind must be either `:rQ` or `:sS`."
    @assert level in [:position, :on_hand] "The policy monitoring level must be either `:position` or `:on_hand`."
    for n in nodes, p in prods
        @assert (n,p) in keys(param1) "The first policy parameter is missing a key for node $n and product $p."
        @assert (n,p) in keys(param2) "The second policy parameter is missing a key for node $n and product $p."
    end

    #initialize action matrix
    action = zeros(length(prods), length(arcs))
    for n in nodes, (k, p) in enumerate(prods)
        if level == :on_hand
            state = filter(i -> i.period == t && i.node == n && i.product == p, env.inv_on_hand).level[1]
        elseif level == :position
            state = filter(i -> i.period == t && i.node == n && i.product == p, env.inv_position).level[1]
        end
        trigger = param1[n,p] > state
        reorder = 0
        if trigger
            if kind == :rQ #rQ policy
                reorder = param2[n,p]
            elseif kind == :sS #sS policy
                reorder = max(param2[n,p] - state, 0)
            end
        else
            reorder = 0
        end

        if supplier_selection == :random
            suppliers = length(inneighbors(env.network, n))
            for src in inneighbors(env.network, n)
                j = findfirst(i -> i == (src, n), arcs)
                action[k, j] = reorder / suppliers #equal split
            end
        elseif supplier_selection == :priority
            supplier_priority = get_prop(env.network, n, :supplier_priority)[p]
            for src in supplier_priority
                if src in env.producers #find available capacity or inventory
                    avail = get_prop(env.network, src, :production_capacity)[p]
                else
                    avail = filter(i -> i.period == env.period && i.node == src && i.product == p, env.inv_on_hand).level[1]
                end
                request = min(avail, reorder) #reorder up to the available amount
                reorder -= request #update reorder quantity for next supplier
                j = findfirst(i -> i == (src, n), arcs) #find index in action matrix
                action[k, j] = request #sate request quantity
            end
        end
    end

    collect(Iterators.flatten(action))
end

function action_space(env::SupplyChainEnv)
    num_products = length(env.products)
    num_edges = ne(env.network)
    num_actions = num_products * num_edges
    ubound = []
    for n in vertices(env.network)
        set_prop!(env.network, n, :max_order, Dict(p => 0. for p in env.products))
    end
    for (source, sink) in Iterators.product(env.producers,env.markets)
        paths = yen_k_shortest_paths(env.network, source, sink, weights(env.network), typemax(Int)).paths
        for path in paths
            capacity = get_prop(env.network, path[1], :production_capacity)
            for i in 2:length(path)
                top_max_order = get_prop(env.network, path[i-1], :max_order)
                if i > 2
                    top_init_inv = get_prop(env.network, path[i-1], :initial_inventory)
                end
                max_order = get_prop(env.network, path[i], :max_order)
                for p in env.products
                    if i == 2
                        max_order[p] += capacity[p]
                    elseif i == 3
                        max_order[p] += top_init_inv[p] + top_max_order[p] * env.num_periods
                    else
                        max_order[p] += top_init_inv[p] + top_max_order[p]
                    end
                end
                set_prop!(env.network, path[i], :max_order, max_order)
            end
        end
    end
    for e in edges(env.network)
        max_order = get_prop(env.network, e.dst, :max_order)
        for p in env.products
            push!(ubound, max_order[p])
        end
    end

    ClosedInterval.(zeros(num_actions), ubound)
end
