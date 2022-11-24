function action_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_edges = ne(env.network)
    return [(0,Inf) for _ in 1:num_products*num_edges]
end

state(env::SupplyChainEnv) = show_state(env).amount

function state_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_nodes = nv(env.network)
    num_arcs = ne(env.network)
    return [(0,Inf) for _ in 1:num_products*(num_nodes*2+num_arcs)] #onhand and unfulfilled for each node and pipeline for each arc
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
    @chain x.inventory begin
        subset( 
            :period => ByRow(==(x.period)),
            :type => ByRow(in(typeset)),
        )
        transform!(
            AsTable([:location,:type]) => ByRow(identity) => :location_type
        )
        unstack(:location_type,:material,:amount,fill=0)
        stack(Not(:location_type), variable_name = :material, value_name = :amount)
        transform!(:location_type => AsTable)
        select!([:type,:location,:material,:amount])
        sort!([:type,:location,:material])
    end
end