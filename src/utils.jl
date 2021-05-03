"""
    check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, bom::Array)

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, bom::Array)
    if bom != [0] #if [0] then single product, no bom
        @assert typeof(bom) <: Matrix{T} where T <: Real "Bill of materials must be a matrix of real numbers."
        @assert size(bom)[1] == size(bom)[2] "Bill of materials must be a square matrix."
    end
    @assert size(bom)[1] == length(mats) "The number of rows and columns in the bill of materials must be equal to the number of materials."
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost]
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :demand_penalty]
    plant_keys = [:production_cost, :production_time, :production_capacity]
    arc_keys = [:sales_price, :transportation_cost, :lead_time]
    for n in nodes, key in all_keys
        @assert key in keys(network.vprops[n]) "$key not stored in distributor node $n."
        for p in mats
            @assert p in keys(network.vprops[n][key]) "Material $p not found in $key on node $n."
            @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
        end
    end
    for n in mrkts, key in market_keys
        @assert key in keys(network.vprops[n]) "$key not stored in market node $n."
        for p in mats
            @assert p in keys(network.vprops[n][key]) "Material $p not found in $key on node $n."
            if key == :demand_frequency
                @assert 0 <= network.vprops[n][key][p] <= 1 "Parameter $key for material $p on node $n must be between 0 and 1."
            elseif key == :demand_distribution
                @assert rand(network.vprops[n][key][p]) isa Number "Parameter $key for material $p on node $n must be a sampleable distribution or an array."
            else
                @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
            end
        end
    end
    for n in plants, key in plant_keys
        @assert key in keys(network.vprops[n]) "$key not stored in producer node $n."
        for p in mats
            @assert p in keys(network.vprops[n][key]) "Material $p not found in $key on node $n."
            @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p on node $n must be non-negative."
        end
    end
    for a in arcs, key in arc_keys
        @assert key in keys(network.eprops[Edge(a...)]) "$key not stored in arc $a."
        if key != :lead_time
            for p in mats
                @assert p in keys(network.eprops[Edge(a...)][key]) "Material $p not found in $key on arc $a."
                @assert network.eprops[Edge(a...)][key][p] >= 0 "Parameter $key for material $p on arc $a must be non-negative."
            end
        else
            @assert rand(network.eprops[Edge(a...)][key]) isa Number "Parameter $key on arc $a must be a sampleable distribution or an array."
        end
    end
    nonsources = [n for n in nodes if !isempty(inneighbors(network, n))]
    for n in nonsources, p in mats
        key = :supplier_priority
        @assert p in keys(network.vprops[n][key]) "Material $p not found in $key on node $n."
        for s in network.vprops[n][key][p]
            @assert s in inneighbors(network, n) "Supplier $s is not a supplier to node $n, but is listed in the supplier priority for that node for material $p."
        end
    end
end
