"""
    reorder_policy(env::SupplyChainEnv, reorder_point::Dict, policy_param::Dict,
                        policy_type::Union{Dict, Symbol} = :rQ, 
                        review_period::Union{Int, AbstractRange, Vector, Dict} = 1,
                        min_order_qty::Union{Real, Dict} = 0,
                        order_multiples::Union{Real, Dict} = 1)

Apply an inventory policy to specify the replinishment orders for each material
    throughout the `SupplyChainEnv`.

# Arguments
- `env`: inventory management environment
- `reorder_point`: the `s` or `r` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_param`: the `S` or `Q` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_type`: `:rQ` for an `(r,Q)` policy, or `:sS` for an `(s,S)` policy. If passing a `Dict`, the policy type should be specified for each node (keys).
- `review_period`: number of periods between each inventory review (Default = `1` for continuous review.). If a `AbstractRange` or `Vector` is used, the `review_period` indicates which periods the review is performed on. If a `Dict` is used, the review period should be specified for each `(node, material)` `Tuple` (keys). The values of this `Dict` can be either `Int`, `AbstractRange`, or `Vector`. Any missing `(node, material)` key will be assigned a default value of 1.
- `min_order_qty`: minimum order quantity (MOQ) at each supply node. If a `Dict` is passed, the MOQ should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of 0.
- `order_multiples`: size increments for each order (default is -1, which means no constraint on order sizes). If a `Dict` is passed, the order multiples should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of -1 (no order multiples enforced).
"""
function reorder_policy(env::SupplyChainEnv, reorder_point::Dict, policy_param::Dict;
                        policy_variable::Union{Dict,Symbol} = :inventory_position,
                        policy_type::Union{Dict, Symbol} = :rQ, 
                        review_period::Union{Int, AbstractRange, Vector, Dict} = 1,
                        min_order_qty::Union{Real, Dict} = 0,
                        order_multiples::Union{Real, Dict} = -1)

    #check review period; if not in review period, send null action
    null_action = zeros(length(env.materials)*ne(env.network))
    !isvalid_period(env.period, review_period) && return null_action

    #read parameters
    nodes = reverse(topological_sort(env.network)) #sort nodes in reverse topological order so that orders are placed moving up the network
    source_nodes = filter(n -> isempty(inneighbors(env.network, n)), nodes) #source nodes (can't place replenishment orders)
    sink_nodes = filter(n -> isempty(outneighbors(env.network, n)), nodes) #sink nodes (can't fulfill replenishment orders)
    request_nodes = setdiff(nodes, source_nodes) #nodes placing requests (all non-source nodes)
    supply_nodes = setdiff(nodes, sink_nodes) #nodes fulfilling requests (all non-sink nodes; used for MOQ)
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    mats = env.materials

    #check inputs
    check_policy_inputs!(reorder_point, policy_param, review_period, min_order_qty, order_multiples, mats, request_nodes, supply_nodes)

    #filter data for policy
    state1_df = filter(:period => j -> j == env.period, env.inventory_position, view=true) #inventory position
    state1_grp = groupby(state1_df, [:node, :material]) #group by node and material
    state2_df = filter(:period => j -> j == env.period, env.echelon_stock, view=true) #inventory position
    state2_grp = groupby(state2_df, [:node, :material]) #group by node and material

    #initialize action matrix
    action = NamedArray(
        zeros(length(mats), length(arcs)),
        (mats, arcs),
        (:material, :arc)
    )
    for n in request_nodes, mat in reverse(get_prop(env.network, n, :node_materials)) #loop materials in reverse topological order (products -> raws)
        #check if an order can be placed
        reorder_point[n,mat] < 0 && continue #if trigger level is negative, skip it
        review_period isa Dict && !isvalid_period(env.period, review_period[n,mat]) && continue #if not a review period, skip
        #review inventory at the node
        state = get_inventory_state(n, mat, policy_variable, state1_grp, state2_grp)
        #if n is a plant, adjust the inventory state of raw materials by that of the expected product orders
        state += get_expected_consumption(env, n, mat, action)
        #check if reorder is triggered & calculate quantity
        if state <= reorder_point[n,mat]
            supplier = get_prop(env.network, n, :supplier_priority)[mat][1] #get first priority supplier
            MOQ = min_order_qty isa Real ? min_order_qty : min_order_qty[(supplier,mat)] #minimum order quantity
            OM = order_multiples isa Real ? order_multiples : order_multiples[(supplier,mat)] #order size
            action[mat, (supplier,n)] = calculate_reorder(n, mat, state, MOQ, OM, policy_type, policy_param)
        end
    end

    return collect(Iterators.flatten(action))
end

"""
    isvalid_period(period::Int, review_period::Union{Int, AbstractRange, Vector, Dict})

Check if current period is a review period.
"""
function isvalid_period(period::Int, review_period::Union{Int, AbstractRange, Vector, Dict})
    if review_period isa Int && !iszero(mod(period,review_period)) 
        return false
    elseif review_period isa Union{AbstractRange,Vector} && !in(period, review_period)
        return false
    else
        return true
    end
end

"""
    check_policy_inputs!(
        reorder_point::Dict, policy_param::Dict, review_period::Union{Int, AbstractRange, Vector, Dict}, 
        min_order_qty::Union{Real, Dict}, order_multiples::Union{Real, Dict}, mats::Vector, request_nodes::Vector, supply_nodes::Vector
    )

Validate inputs to inventory policy.
"""
function check_policy_inputs!(
    reorder_point::Dict, policy_param::Dict, review_period::Union{Int, AbstractRange, Vector, Dict}, 
    min_order_qty::Union{Real, Dict}, order_multiples::Union{Real, Dict}, 
    mats::Vector, request_nodes::Vector, supply_nodes::Vector
)
    #check inputs
    for mat in mats
        for n in request_nodes
            for param in [reorder_point, policy_param] 
                if !in((n,mat), keys(param)) 
                    param[(n,mat)] = -1 #if no policy given for a node/material, set the params to -1 so that a reorder is never triggered
                end
            end
            if review_period isa Dict && !in((n,mat), keys(review_period))
                review_period[(n,mat)] = 1 #set to default value of 1
            end
        end
        for n in supply_nodes
            if min_order_qty isa Dict && !in((n,mat), keys(min_order_qty))
                min_order_qty[(n,mat)] = 0 #set to default value of 0
            end
            if order_multiples isa Dict && !in((n,mat), keys(order_multiples))
                order_multiples[(n,mat)] = -1 #set to default value of -1
            end
        end
    end
end

"""
    get_inventory_state(
        n::Int, mat::Union{Symbol,String}, policy_variable::Union{Dict,Symbol},
        state1_grp::GroupedDataFrame, state2_grp::GroupedDataFrame
    )

Check the inventory state at node `n` for material `mat`.
"""
function get_inventory_state(
    n::Int, mat::Union{Symbol,String}, policy_variable::Union{Dict,Symbol},
    state1_grp::GroupedDataFrame, state2_grp::GroupedDataFrame
)
    pol_var = policy_variable isa Dict ? policy_variable[n] : policy_variable
    @assert pol_var in [:inventory_position, :echelon_stock] "Policy variable must be either `:inventory_position` or `:echelon_stock`."
    if pol_var == :inventory_position
        state = state1_grp[(node = n, material = mat)].level[1]
    elseif pol_var == :echelon_stock
        state = state2_grp[(node = n, material = mat)].level[1]
    end

    return state
end

"""
    get_expected_consumption(env::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, action::NamedArray)

Get expected material consumption at plants.
"""
function get_expected_consumption(env::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, action::NamedArray)
    consume = 0
    if isconsumed(env, n, mat)
        bom = get_prop(env.network, n, :bill_of_materials)
        mprods = names(filter(i -> i < 0, bom[mat,:]), 1) #find which products are produced from mat
        for mprod in mprods
            ordered = action[mprod, (n, n)] #amount of product ordered
            stoich = bom[mat, mprod] #stoichiometry
            consume += ordered * stoich #amount that will be consumed to fulfill the order (negative)
        end
    end

    return consume
end

"""
    calculate_reorder(
        n::Int, mat::Union{Symbol,String}, state::Float64, MOQ::Real, OM::Real,
        policy_type::Union{Dict, Symbol}, policy_param::Dict
    )

Calculate reorder quantity for material `mat` at node `n`.
"""
function calculate_reorder(
    n::Int, mat::Union{Symbol,String}, state::Float64, MOQ::Real, OM::Real,
    policy_type::Union{Dict, Symbol}, policy_param::Dict
)
    pol_type = policy_type isa Dict ? policy_type[n] : policy_type #policy type
    @assert pol_type in [:rQ, :sS] "The policy type must be either `:rQ` or `:sS`."

    #get order size
    if pol_type == :rQ #rQ policy
        Q = policy_param[n,mat] #order size
    elseif pol_type == :sS #sS policy
        Q = max(policy_param[n,mat] - state, 0) #order size (only order if state is below S)
    end

    #apply MOQ and OM constraints (only if order is positive)
    if Q > 0 
        if OM < 0 #no constraint on order multiples
            reorder = max(Q, MOQ)
        else #apply order size constraint
            reorder = MOQ + OM*ceil((Q-MOQ)/OM) 
        end
    else
        reorder = 0
    end
    
    return reorder
end

"""
    simulate_policy!(env::SupplyChainEnv, args...)

Step through a simulation using a specified reorder policy. `args` are the
arguments that are passed to the `reorder_policy` function.
"""
function simulate_policy!(env::SupplyChainEnv, args...; kwargs...)
    for _ in 1:env.num_periods
        action = reorder_policy(env, args...; kwargs...)
        (env)(action)
    end

    calculate_service_measures!(env)
end
