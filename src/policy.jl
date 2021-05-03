"""
    reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                level::Symbol = :position, kind::Symbol = :rQ,
                supplier_selection::Symbol = :priority, review_period::Int = 1)

Apply an inventory policy to specify the replinishment orders for each material
    throughout the `SupplyChainEnv`.
"""
function reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                level::Symbol = :position, kind::Symbol = :rQ,
                supplier_selection::Symbol = :priority, review_period::Int = 1)

    #read parameters
    t = env.period
    nodes = [n for n in vertices(env.network) if !isempty(inneighbors(env.network, n))] #all non-source nodes can place orders
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    mats = env.materials

    #check inputs
    @assert kind in [:rQ, :sS] "The policy kind must be either `:rQ` or `:sS`."
    @assert level in [:position, :on_hand] "The policy monitoring level must be either `:position` or `:on_hand`."
    for n in nodes, p in mats
        @assert (n,p) in keys(param1) "The first policy parameter is missing a key for node $n and material $p."
        @assert (n,p) in keys(param2) "The second policy parameter is missing a key for node $n and material $p."
    end

    #check review period
    if !iszero(mod(t,review_period)) #if not in review period, send null action
        return zeros(length(mats)*length(arcs))
    end

    #initialize action matrix
    action = zeros(length(mats), length(arcs))
    for n in nodes, (k, p) in enumerate(mats)
        param1[n,p] < 0 && continue #if trigger level is negative, skip it
        trigger = false #trigger an inventory replenishment order
        reorder = 0
        #check if reorder is triggered
        if level == :on_hand
            state = filter(i -> i.period == t && i.node == n && i.material == p, env.inv_on_hand).level[1]
            if t-review_period <= 0
                trigger = param1[n,p] >= state #trigger if initial inventory is below trigger point
            else
                state0 = filter(i -> i.period == t-review_period && i.node == n && i.material == p, env.inv_on_hand).level[1]
                trigger = state0 > param1[n,p] >= state #only reorder first time going past the threshold
            end
        elseif level == :position
            state = filter(i -> i.period == t && i.node == n && i.material == p, env.inv_position).level[1]
            trigger = param1[n,p] >= state
        end
        #set reorder policy
        if trigger
            if kind == :rQ #rQ policy
                reorder = param2[n,p]
            elseif kind == :sS #sS policy
                reorder = max(param2[n,p] - state, 0)
            end
        else
            continue
        end
        #assign reorder quantities
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
                    avail = filter(i -> i.period == env.period && i.node == src && i.material == p, env.inv_on_hand).level[1]
                end
                request = min(avail, reorder) #reorder up to the available amount
                reorder -= request #update reorder quantity for next supplier
                j = findfirst(i -> i == (src, n), arcs) #find index in action matrix
                action[k, j] = request #sate request quantity
            end
        end
    end

    return collect(Iterators.flatten(action))
end

"""
    simulate_policy!(env::SupplyChainEnv, args...)

Step through a simulation using a specified reorder policy. `args` are the
arguments that are passed to the `reorder_policy` function.
"""
function simulate_policy!(env::SupplyChainEnv, args...)
    for t in 1:env.num_periods
        action = reorder_policy(env, args...)
        (env)(action)
    end
end