"""
    reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,
                        policy_type::Union{Dict, Symbol} = :rQ, 
                        review_period::Union{Int, StepRange, Vector, Dict} = 1,
                        min_order_qty::Union{Real, Dict} = 0)

Apply an inventory policy to specify the replinishment orders for each material
    throughout the `SupplyChainEnv`.

# Arguments
- `env::SupplyChainEnv`: inventory management environment
- `param1::Dict`: the `s` or `r` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `param2::Dict`: the `S` or `Q` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_type::Union{Symbol, Dict}`: `:rQ` for an `(r,Q)` policy, or `:sS` for an `(s,S)` policy. If passing a `Dict`, the policy type should be specified for each node (keys).
- `review_period::Union{Int, StepRange, Vector, Dict}`: number of periods between each inventory review (Default = `1` for continuous review.). If a `StepRange` or `Vector` is used, the `review_period` indicates which periods the review is performed on. If a `Dict` is used, the review period should be specified for each `(node, material)` `Tuple` (keys). The values of this `Dict` can be either `Int`, `StepRange`, or `Vector`. Any missing `(node, material)` key will be assigned a default value of 1.
- `min_order_qty::Union{Real, Dict}`: minimum order quantity (MOQ) at each supply node. If a `Dict` is passed, the MOQ should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of 0.
"""
function reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict;
                        policy_variable::Union{Dict,Symbol} = :inventory_position,
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
    state1_df = filter(:period => j -> j == env.period, env.inventory_position, view=true) #inventory position
    state1_grp = groupby(state1_df, [:node, :material]) #group by node and material
    state2_df = filter(:period => j -> j == env.period, env.echelon_stock, view=true) #inventory position
    state2_grp = groupby(state2_df, [:node, :material]) #group by node and material

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
        pol_var = policy_variable isa Dict ? policy_variable[n] : policy_variable
        if pol_var == :inventory_position
            state = state1_grp[(node = n, material = p)].level[1]
        elseif pol_var == :echelon_stock
            state = state2_grp[(node = n, material = p)].level[1]
        else
            @assert pol_var in [:inventory_position, :echelon_stock] "Policy variable must be either `:inventory_position` or `:echelon_stock`."
        end
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
function simulate_policy!(env::SupplyChainEnv, args...; kwargs...)
    for t in 1:env.num_periods
        action = reorder_policy(env, args...; kwargs...)
        (env)(action)
    end
    
    calculate_service_measures!(env)
end
