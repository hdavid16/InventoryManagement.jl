"""
    show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)

Visualize a replenishment order vector as a DataFrame indicating how much of
which material is being requested on each arc.
"""
function show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)
    mats = env.materials
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    act = reshape(action, (length(mats), length(arcs)))
    return DataFrame(:material => mats, 
                     [Symbol(a) => act[:,i] for (i,a) in enumerate(arcs)]...)
end

"""
    check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, bom::Array,
                    num_periods::Int)

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, bom::Array,
                    num_periods::Int)
    if bom != [0] #if [0] then single product, no bom
        @assert typeof(bom) <: Matrix{T} where T <: Real "Bill of materials must be a matrix of real numbers."
        @assert size(bom)[1] == size(bom)[2] "Bill of materials must be a square matrix."
    end
    @assert size(bom)[1] == length(mats) "The number of rows and columns in the bill of materials must be equal to the number of materials."
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost]
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :demand_penalty, :demand_sequence]
    plant_keys = [:production_cost, :production_time, :production_capacity]
    arc_keys = [:sales_price, :transportation_cost, :lead_time]
    for n in nodes, key in all_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set its value to 0
                network.vprops[n][key][p] = 0.
            end
            @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
        end
    end
    for n in mrkts, key in market_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set its value to 0
                if key == :demand_distribution
                    tmp = Dict(p => [0])
                    network.vprops[n][key] = merge(network.vprops[n][key], tmp)
                elseif key == :demand_sequence
                    network.vprops[n][key][p] = zeros(num_periods)
                else
                    network.vprops[n][key][p] = 0.
                end
            end
            if key == :demand_frequency
                @assert 0 <= network.vprops[n][key][p] <= 1 "Parameter $key for material $p on node $n must be between 0 and 1."
            elseif key == :demand_distribution
                @assert rand(network.vprops[n][key][p]) isa Number "Parameter $key for material $p on node $n must be a sampleable distribution or an array."
            elseif key == :demand_sequence
                @assert length(network.vprops[n][key][p]) == num_periods "The demand sequence for material $p at node $n must be a vector with $num_periods entries."
            else
                @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
            end
        end
    end
    for n in plants, key in plant_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set its value to 0
                network.vprops[n][key][p] = 0.
            end
            @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
        end
    end
    for a in arcs, key in arc_keys
        !in(key, keys(network.eprops[Edge(a...)])) && set_prop!(network, Edge(a...), key, Dict()) #create empty params for arcs if not specified
        if key != :lead_time
            for p in mats
                if !in(p, keys(network.eprops[Edge(a...)][key])) #if material not specified, add it to the dict and set its value to 0
                    network.eprops[Edge(a...)][key][p] = 0.
                end
                @assert network.eprops[Edge(a...)][key][p] >= 0 "Parameter $key for material $p on arc $a must be non-negative."
            end
        else
            @assert rand(network.eprops[Edge(a...)][key]) isa Number "Parameter $key on arc $a must be a sampleable distribution or an array."
        end
    end
    nonsources = [n for n in nodes if !isempty(inneighbors(network, n))]
    for n in nonsources, p in mats
        key = :supplier_priority
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
        if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and add the node's predecessors
            network.vprops[n][key][p] = inneighbors(network, n)
        end
        for s in network.vprops[n][key][p]
            @assert s in inneighbors(network, n) "Supplier $s is not a supplier to node $n, but is listed in the supplier priority for that node for material $p."
        end
    end
end
