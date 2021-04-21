function (x::SupplyChainEnv)(action)
    #validate action input (i.e., type = DataFrame, fields = ["product" => String, "Arc_i_j (same order as edges(network))" => Float64])

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

    #place requests
    for i in 1:nrow(action)
        p = action[i, :product] #product requested
        for a in edges(x.network)
            amount = copy(action[i, Symbol((a.src,a.dst))]) #amount requested
            prod_time = 0 #initialize production time
            #accept or adjust requests
            if a.src in x.producers
                capacity = get_prop(x.network, a.src, :production_capacity)[p] #get production capacity
                accepted = min(amount, capacity) #accepted request
                action[i, Symbol((a.src,a.dst))] = accepted #overwrite request to comply with constraints
                if amount > capacity
                    @warn "Replenishment request for product $p to node $(a.src) was reduced by $(amount - action[i,Symbol((a.src,a.dst))]) due to insufficient production capacity."
                end
                if accepted > 0 #schedule production
                    prod_time = get_prop(x.network, a.src, :production_time)[p] #production time
                    push!(x.production, ((a.src,a.dst), p, accepted, prod_time))
                end
            else
                supply = filter(i -> i.period == x.period && i.node == a.src && i.product == p, x.inv_on_hand).level[1] #on_hand inventory at supplier
                accepted = min(amount, supply) #accepted request
                action[i, Symbol((a.src,a.dst))] = accepted #overwrite request to comply with constraints
                if amount > supply
                    @warn "Replenishment request for product $p to node $(a.src) was reduced by $(amount - action[i,Symbol((a.src,a.dst))]) due to insufficient inventory."
                end
                if accepted > 0 #subtract sent onhand inventory, add sent pipeline inventory
                    x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                                  (x.inv_on_hand.node .== a.src) .&
                                  (x.inv_on_hand.product .== p), :level] .-= accepted
                    x.inv_pipeline[(x.inv_pipeline.period .== x.period) .&
                                   (string.(x.inv_pipeline.arc) .== string(a)) .&
                                   (x.inv_pipeline.product .== p), :level] .+= accepted
                end
            end
            #store shipments and lead times
            if accepted > 0
                lead = leads[a] + prod_time
                push!(x.shipments, [(a.src,a.dst), p, accepted, lead]) #update active shipments
            else
                lead = 0 #no request made
            end
            push!(x.replenishments, [x.period, (a.src,a.dst), p, accepted, lead])
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

function policy(env::SupplyChainEnv, param::DataFrame, level::Symbol = :position)
    t = env.period
    nodes = union(env.distributors, env.markets)
    action = DataFrame(:product => env.products)
    for n in nodes
        if level == :on_hand
            state = filter(i -> i.period == t && i.node == n, env.inv_on_hand)[:,["product", "level"]]
        elseif level == :position
            state = filter(i -> i.period == t && i.node == n, env.inv_position)[:,["product", "level"]]
        end
        df = leftjoin(param, state, on="product")
        df.trigger = df[:,2] .> df.level #detect trigger on inventory
        ub = names(param)[3]
        df.reorder = zeros(nrow(df))
        if ub == "Q" #rQ policy
            df.reorder = df.trigger .* df.Q
        elseif ub == "S" #sS policy
            df.reorder = df.trigger .* (df.S .- df.level)
        end

        suppliers = length(inneighbors(env.network, n))
        for src in inneighbors(env.network, n)
            action[:, Symbol((src,n))] = df.reorder / suppliers #equal split
        end
    end

    return action
end
