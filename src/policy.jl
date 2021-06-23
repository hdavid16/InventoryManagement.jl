"""
    reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                        kind::Symbol = :rQ, review_period::Int = 1)

Apply an inventory policy to specify the replinishment orders for each material
    throughout the `SupplyChainEnv`.
"""
function reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                        kind::Symbol = :rQ, review_period::Int = 1)

    #check review period
    if !iszero(mod(env.period,review_period)) #if not in review period, send null action
        return zeros(length(env.materials)*ne(env.network))
    end

    #read parameters
    t = env.period
    nodes = [n for n in vertices(env.network) if !isempty(inneighbors(env.network, n))] #all non-source nodes can place orders
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    mats = env.materials

    #check inputs
    @assert kind in [:rQ, :sS] "The policy kind must be either `:rQ` or `:sS`."
    for n in nodes, p in mats, param in [param1, param2] #if no policy given for a node/material, set the params to -1 so that a reorder is never triggered
        if !in((n,p), keys(param))
            param[(n,p)] = -1
        end
    end

    #initialize action matrix
    action = zeros(length(mats), length(arcs))
    state_df = filter([:period] => j -> j == t, env.inv_position, view=true) #on_hand inventory
    state_grp = groupby(state_df, [:node, :material]) #group by node and material
    for n in nodes, (k, p) in enumerate(mats)
        param1[n,p] < 0 && continue #if trigger level is negative, skip it
        reorder = 0
        #check if reorder is triggered & set reorder policy
        state = state_grp[(node = n, material = p)].level[1]
        if state <= param1[n,p]
            if kind == :rQ #rQ policy
                reorder = param2[n,p]
            elseif kind == :sS #sS policy
                reorder = max(param2[n,p] - state, 0)
            end
        else
            continue
        end
        #assign reorder quantity
        supplier = get_prop(env.network, n, :supplier_priority)[p][1]
        j = findfirst(i -> i == (supplier, n), arcs) #find index in action matrix
        action[k, j] = reorder #sate request quantity
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
        mod(env.period,24*7) == 0 && println("Week $(env.period/24/7)")
    end
end
