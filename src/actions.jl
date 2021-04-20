function (x::SupplyChainEnv)(action)
    #validate action input (i.e., type = DataFrame, fields = ["product" => String, "Arc_i_j (same order as edges(network))" => Float64])

    #increase period counter
    x.period += 1
    #move active shipments forward one period
    x.shipments.lead .-= 1
    #sample lead times
    leads = Dict(a => rand(get_prop(x.network, a, :lead_time)) for a in edges(x.network)) #sample lead time distribution
    #place requests
    for i in 1:nrow(action)
        p = action[i, :product] #product requested NOTE: must be sorted same as params
        for a in edges(x.network)
            amount = copy(action[i, Symbol((a.src,a.dst))]) #amount requested
            #accept or adjust requests
            if a.src in x.producers
                capacity = filter(j -> j.product == p, get_prop(x.network, a.src, :params)).production_capacity[1]
                action[i, Symbol((a.src,a.dst))] = min(amount, capacity)
                if amount > capacity
                    @warn "Replenishment request for product $p to node $(a.src) was reduced by $(amount - action[i,Symbol((a.src,a.dst))]) due to insufficient production capacity."
                end
            else
                supply = filter(i -> i.period == x.period && i.node == a.src && i.product == p, x.inv_on_hand).level[1]
                action[i, Symbol((a.src,a.dst))] = min(amount, supply)
                if amount > supply
                    @warn "Replenishment request for product $p to node $(a.src) was reduced by $(amount - action[i,Symbol((a.src,a.dst))]) due to insufficient inventory."
                end
            end
            #store shipments and lead times
            if action[i, Symbol((a.src,a.dst))] > 0
                lead = leads[a]
                push!(x.shipments, [(a.src,a.dst), p, action[i, Symbol((a.src,a.dst))], lead]) #update active shipments
            else
                lead = 0 #no request made
            end
            push!(x.replenishments, [x.period, (a.src,a.dst), p, action[i, Symbol((a.src,a.dst))], lead])
        end
    end
    #update non market inventories
    arrivals = filter(i -> iszero(i.lead), x.shipments) #find active shipments with 0 lead time
    departures = filter(i -> i.period == x.period, x.replenishments) #requests sent in current period
    for i in 1:nrow(departures)
        arc = departures[i,"arc"] #arc sending inventory
        node = arc[1] #(source) node sending inventory
        p = departures[i,"product"] #product sent

        #update on-hand inventory
        sent = departures[i,"amount"] #amount sent
        received = filter(j -> j.arc[end] == node && j.product == p, arrivals).amount #amount received from incoming shipment
        if isempty(received) #if none is received, set to 0
            received = [0]
        end
        prev = filter(j -> j.period == x.period-1 && j.node == node && j.product == p, x.inv_on_hand).level[1] #previous inventory level
        if node in x.producers
            onhand = 0
        else
            onhand = prev + received[1] - sent #new on-hand inventory
        end
        push!(x.inv_on_hand, [x.period, node, p, onhand]) #update inventory

        #update pipeline inventory
        leaving = filter(j -> j.arc == arc && j.product == p, arrivals).amount #amount leaving the pipeline
        if isempty(leaving) #if non leaves pipeline, set to 0
            leaving = [0]
        end
        prev = filter(j -> j.period == x.period-1 && j.arc == arc && j.product == p, x.inv_pipeline).level[1] #previous pipeline inv level
        push!(x.inv_pipeline, [x.period, arc, p, prev + sent - leaving[1]]) #update inventory

        #update inventory position
        upstream = sum(filter(j -> j.period == x.period && j.arc[end] == node && j.product == p, x.inv_pipeline).level)
        push!(x.inv_position, [x.period, node, p, onhand + upstream]) #update inventory
    end
    #update market inventories
    for node in x.markets, p in x.products
        prev = filter(j -> j.period == x.period-1 && j.node == node && j.product == p, x.inv_on_hand).level[1] #previous inventory level
        push!(x.inv_on_hand, [x.period, node, p, prev]) #intialize with previous inventory levels
    end
    arrivals = filter(i -> iszero(i.lead) && i.arc[end] in x.markets, x.shipments) #shipments arriving at sink nodes (markets)
    for i in 1:nrow(arrivals)
        node = arrivals[i,:arc][end] #market node
        p, amount = arrivals[i,["product","amount"]]
        x.inv_on_hand[(x.inv_on_hand.period .== x.period) .&
                      (x.inv_on_hand.node .== node) .&
                      (x.inv_on_hand.product .== p), :level] .+= amount #update inventory
    end
    filter!(i -> i.lead > 0, x.shipments) #remove shipments that arrived

    #markets open and demand occurs
    for m in x.markets
        freq = Bernoulli.(get_prop(x.network, m, :params).demand_frequency)
        dmnd = get_prop(x.network, m, :params).market_demand
        quant = rand.(values(freq)) .* rand.(values(dmnd))
        for (q, p) in zip(quant, x.products)
            demand = [x.period, m, p, q, 0, 0, 0]
            inv = filter(i -> i.period == x.period && i.node == m && i.product == p, x.inv_on_hand).level[1]
            if x.backlog
                q += filter(i -> i.period == x.period-1 && i.market == m && i.product == p, x.demand).backlog[1]
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
    end

    #calculate profit at each node
    for n in vertices(x.network)
        profit = 0 #initialize node profit
        #pay inventory holding cost
        params0 = get_prop(x.network, n, :params)[:,["product", "holding_cost", "production_cost", "market_price", "demand_penalty"]]
        onhand = filter(j -> j.period == x.period && j.node == n, x.inv_on_hand)[:,["product", "level"]]
        df = leftjoin(params0, onhand, on="product")
        profit -= sum(skipmissing(df.level .* df.holding_cost))
        #pay suppliers for received inventory and pay transportation cost
        for pred in inneighbors(x.network, n)
            params = get_prop(x.network, pred, n, :params)[:,["product", "sales_price", "transportation_cost"]]
            #pay purchase of inventory
            purchased = filter(j -> j.arc == (pred, n), arrivals)[:,["product", "amount"]]
            if !isempty(purchased)
                df = leftjoin(params, purchased, on="product")
                profit -= sum(skipmissing(df.amount .* df.sales_price))
            end
            #pay transportation (pipeline hoding) cost
            intransit = filter(j -> j.period == x.period && j.arc == (pred, n), x.inv_pipeline)[:,["product", "level"]]
            df = leftjoin(params, intransit, on="product")
            profit -= sum(skipmissing(df.level .* df.transportation_cost))
        end
        #receive payment for delivered inventory and pay production costs
        for succ in outneighbors(x.network, n)
            params = get_prop(x.network, n, succ, :params)[:,["product", "sales_price"]]
            #receive payment for delivered inventory
            sold = filter(j -> j.arc == (n, succ), arrivals)[:,["product", "amount"]]
            if !isempty(sold)
                df = leftjoin(params, sold, on="product")
                profit += sum(skipmissing(df.amount .* df.sales_price))
            end
            #pay production cost
            if n in x.producers
                produced = filter(j -> j.period == x.period && j.arc == (n, succ), x.replenishments)
                df1 = leftjoin(params0, produced, on="product")
                profit -= sum(skipmissing(df1.amount .* df1.production_cost))
            end
        end
        #sales profit at end distributors (and penalize for unfulfilled demand)
        if n in x.markets
            sold = filter(j -> j.period == x.period && j.market == n, x.demand)
            df1 = leftjoin(params0, sold, on="product")
            profit += sum(skipmissing(df1.sale .* df1.market_price))
            profit -= sum(skipmissing(df1.unfulfilled .* df1.demand_penalty))
        end
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
