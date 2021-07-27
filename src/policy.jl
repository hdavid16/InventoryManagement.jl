"""
    reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                        policy_type::Union{Dict, Symbol} = :rQ, 
                        review_period::Union{Int, StepRange, Vector, Dict} = 1,
                        min_order_qty::Union{Int, Dict} = 0)

Apply an inventory policy to specify the replinishment orders for each material
    throughout the `SupplyChainEnv`.

If `review_period` or `min_order_qty` are `Dict`s, they must have (node, material) `Tuples` as the keys.
"""
function reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                        policy_type::Union{Dict, Symbol} = :rQ, 
                        review_period::Union{Int, StepRange, Vector, Dict} = 1,
                        min_order_qty::Union{Real, Dict} = 0)

    #check review period
    null_action = zeros(length(env.materials)*ne(env.network))
    if review_period isa Int && !iszero(mod(env.period,review_period)) #if not in review period, send null action
        return null_action
    elseif review_period isa Union{StepRange,Vector} && !in(env.period, review_period)
        return null_action
    end

    #read parameters
    nodes = [n for n in vertices(env.network) if !isempty(inneighbors(env.network, n))] #all non-source nodes can place orders
    supply_nodes = [n for n in vertices(env.network) if !isempty(outneighbors(env.network, n))] #all non-sink nodes can place orders
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    mats = env.materials

    #check inputs
    for p in mats
        for n in nodes
            for param in [param1, param2] 
                if !in((n,p), keys(param)) 
                    param[(n,p)] = -1 #if no policy given for a node/material, set the params to -1 so that a reorder is never triggered
                end
            end
            if review_period isa Dict && !in((n,p), keys(review_period))
                review_period[(n,p)] = 1 #set to default value of 1
            end
        end
        for n in supply_nodes
            if min_order_qty isa Dict && !in((n,p), keys(min_order_qty))
                min_order_qty[(n,p)] = 0 #set to default value of 0
            end
        end
    end

    #filter data for policy
    state_df = filter(:period => j -> j == env.period, env.inv_position, view=true) #inventory position
    state_grp = groupby(state_df, [:node, :material]) #group by node and material

    #create action matrix
    action = zeros(length(mats), length(arcs)) #initialize action matrix
    for n in nodes, (k, p) in enumerate(mats)
        #if trigger level is negative, skip it
        param1[n,p] < 0 && continue 
        #if not a review period, skip
        if review_period isa Dict
            review_period[n,p] isa Int && !iszero(mod(env.period, review_period[n,p])) && continue #trigger if not in review period
            review_period[n,p] isa Union{StepRange,Vector} && !in(env.period, review_period[n,p]) && continue
        end
        #get MOQ
        supplier = get_prop(env.network, n, :supplier_priority)[p][1] #get first priority supplier
        MOQ = min_order_qty isa Real ?
                min_order_qty :
                min_order_qty[(supplier,p)]
        #initialize reorder quantity
        reorder = 0
        #check if reorder is triggered & set reorder policy
        state = state_grp[(node = n, material = p)].level[1]
        if state <= param1[n,p]
            pol_type = policy_type isa Dict ? policy_type[n] : policy_type
            if pol_type == :rQ #rQ policy
                reorder = max(param2[n,p], MOQ)
            elseif pol_type == :sS #sS policy
                reorder = state < param2[n,p] ? max(param2[n,p] - state, MOQ) : 0
            else
                @assert pol_type in [:rQ, :sS] "The policy type must be either `:rQ` or `:sS`."
            end
        else
            continue
        end
        #assign reorder quantity
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
