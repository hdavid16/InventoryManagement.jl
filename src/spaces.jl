function action_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_edges = ne(env.network)
    return [Interval{:closed,:open}(0,Inf) for _ in num_products*num_edges]
end

function state(env::SupplyChainEnv)
    typeset = Set([:on_hand,:pipeline,:unfulfilled])
    return sort(
        subset(env.inventory, 
            :period => ByRow(==(env.period)),
            :type => ByRow(in(typeset)),
            view=true
        ),
        [:type,:location,:material]
    ).amount
end

function state_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_nodes = nv(env.network)
    num_types = 3 #onhand, pipeline, unfulfilled
    return [Interval{:closed,:open}(0,Inf) for _ in num_products*num_nodes*num_types]
end

"""
    show_action(x::SupplyChainEnv, action::Vector{T} where T <: Real)

Convert a replenishment order vector into a NamedArray indicating how much of
each material (rows) is being requested on each arc (columns).
"""
function show_action(x::SupplyChainEnv, action::Vector{T} where T <: Real)
    mats = x.materials
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    return NamedArray(
        reshape(action, (length(mats), length(arcs))),
        (x.materials, arcs),
        (:material, :arc)
    )
end

"""
    show_state(x::SupplyChainEnv)

Show the system state as a DataFrame.
"""
function show_state(x::SupplyChainEnv)
    typeset = Set([:on_hand,:pipeline,:unfulfilled])
    return sort(
        subset(x.inventory, 
            :period => ByRow(==(x.period)),
            :type => ByRow(in(typeset)),
            view=true
        ),
        [:type,:location,:material]
    )
end