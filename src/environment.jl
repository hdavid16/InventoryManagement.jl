abstract type AbstractEnv end

"""
    Supply Chain Simulation Environment Constructor

# Fields:
- `network::MetaDiGraph`: Supply chain network directed graph.
- `markets::Vector`: Vector of market nodes where demand occurs.
- `producers::Vector`: Vector of producer nodes where material transformation occurs.
- `distributors::Vector`: Vector of distribution centers (excludes market nodes).
- `materials::Vector`: Vector with the names of all materials in the system.
- `bill_of_materials::Array`: Square matrix with BOM (rows = inputs, cols = ouputs); indices follow materials list; row with positive value is a co-product, row with negative value is a input (reactant).
- `inv_on_hand::DataFrame`: Timeseries with on hand inventories @ each node.
- `inv_pipeline::DataFrame`: Timeseries with pipeline inventories on each arc.
- `inv_position::DataFrame`: Timeseries with inventory positions @ each node.
- `replenishments::DataFrame`: Timeseries with replenishment orders placed on each arc.
- `shipments::DataFrame`: Temp table with active shipments and time to arrival on each arc.
- `production::DataFrame`: Temp table with active material production commited to an arc and time to ship.
- `demand::DataFrame`: Timeseries with realization of demand at each market, and amounts sold, unfulfilled demand, and backlog.
- `profit::DataFrame`: Timeseries with profit @ each node.
- `reward::Float64`: Final reward in the system (used for RL)
- `period::Int`: Current period in the simulation.
- `num_periods::Int`: Number of periods in the simulation.
- `discount::Float64`: Time discount factor (i.e. interest rate).
- `backlog::Bool`: Indicator if backlogging is allowed.
- `reallocate::Bool`: Indicator if unfulfilled requests should be reallocated to alternate suppliers.
- `seed::Int`: Random seed.
"""
mutable struct SupplyChainEnv <: AbstractEnv
    network::MetaDiGraph
    markets::Vector
    producers::Vector
    distributors::Vector
    materials::Vector
    bill_of_materials::Array
    inv_on_hand::DataFrame
    inv_pipeline::DataFrame
    inv_position::DataFrame
    replenishments::DataFrame
    shipments::DataFrame
    production::DataFrame
    demand::DataFrame
    profit::DataFrame
    reward::Float64
    period::Int
    num_periods::Int
    discount::Float64
    backlog::Bool
    reallocate::Bool
    seed::Int
end

"""
    SupplyChainEnv(network::MetaDiGraph, num_periods::Int;
                        discount::Float64=0.0, backlog::Bool=true,
                        reallocate::Bool=true,
                        seed::Int=0)

Create a `SupplyChainEnv` from a directed graph with metadata (`MetaDiGraph`).
"""
function SupplyChainEnv(network::MetaDiGraph, num_periods::Int;
                        discount::Float64=0.0, backlog::Bool=true,
                        reallocate::Bool=true,
                        seed::Int=0)
    #get main nodes
    nodes = vertices(network)
    #get main edges
    arcs = [(e.src,e.dst) for e in edges(network)]
    #get end distributors, producers, and distribution centers
    mrkts = [n for n in nodes if isempty(outneighbors(network, n))] #markets must be sink nodes
    plants = [n for n in nodes if :production_cost in keys(network.vprops[n])]
    dcs = setdiff(nodes, mrkts, plants)
    #get materials
    mats = get_prop(network, :materials)
    bom = get_prop(network, :bill_of_materials)
    #check inputs
    check_inputs(network, nodes, arcs, mrkts, plants, mats, bom, num_periods)
    #create logging dataframes
    inv_on_hand = DataFrame(:period => Int[], :node => Int[], :material => [], :level => Float64[], :discarded => [])
    for n in nodes, p in mats
        init_inv = get_prop(network, n, :initial_inventory)
        push!(inv_on_hand, (0, n, p, init_inv[p], 0))
    end
    inv_pipeline = DataFrame(:period => zeros(Int, length(mats)*length(arcs)),
                             :arc => repeat(arcs, inner = length(mats)),
                             :material => repeat(mats, outer = length(arcs)),
                             :level => zeros(length(mats)*length(arcs)))
    inv_position = select(inv_on_hand, [:period, :node, :material, :level])
    replenishments = DataFrame(:period => Int[],#zeros(Int, length(mats)*length(arcs)),
                               :arc => Tuple[],#repeat(arcs, inner = length(mats)),
                               :material => Any[],#repeat(mats, outer = length(arcs)),
                               :amount => Float64[],#zeros(length(mats)*length(arcs)),
                               :lead => Int[],#zeros(Int, length(mats)*length(arcs)),
                               :unfulfilled => Float64[],#zeros(length(mats)*length(arcs)),
                               :reallocated => Any[])#Any[missing for i in 1:length(mats)*length(arcs)])
    shipments = DataFrame(:arc => [],
                          :material => [],
                          :amount => Float64[],
                          :lead => Int[])
    production = DataFrame(:arc => [],
                           :material => [],
                           :amount => Float64[],
                           :lead => Int[])
    demand = DataFrame(:period => Int[],#zeros(Int, length(mats)*length(mrkts)),
                       :node => Int[],#repeat(mrkts, inner = length(mats)),
                       :material => Any[],#repeat(mats, outer = length(mrkts)),
                       :demand => Float64[],#zeros(length(mats)*length(mrkts)),
                       :sale => Float64[],#zeros(length(mats)*length(mrkts)),
                       :unfulfilled => Float64[])#zeros(length(mats)*length(mrkts)))
    profit = DataFrame(:period => Int[],#zeros(Int, length(nodes)),
                       :value => Float64[],#zeros(length(nodes)),
                       :node => Int[])#nodes)
    reward = 0
    period = 0
    num_periods = num_periods
    env = SupplyChainEnv(network, mrkts, plants, dcs, mats, bom, inv_on_hand, inv_pipeline, inv_position,
                    replenishments, shipments, production, demand,
                    profit, reward, period, num_periods, discount, backlog, reallocate, seed)

    Random.seed!(env)

    return env
end

"""
    reset!(env::SupplyChainEnv)

Reset a `SupplyChainEnv` (empty all logging dataframes and set simulation time to 0).
"""
function reset!(env::SupplyChainEnv)
    env.period = 0
    env.reward = 0
    filter!(i -> i.period == 0, env.inv_on_hand)
    filter!(i -> i.period == 0, env.inv_pipeline)
    filter!(i -> i.period == 0, env.inv_position)
    filter!(i -> i.period == 0, env.replenishments)
    filter!(i -> i.period == 0, env.demand)
    filter!(i -> i.period == 0, env.profit)
    env.shipments = DataFrame(:arc => [],
                          :material => [],
                          :amount => Float64[],
                          :lead => Int[])
    env.production = DataFrame(:arc => [],
                           :material => [],
                           :amount => Float64[],
                           :lead => Int[])

    Random.seed!(env)
end

"""
    is_terminated(env::SupplyChainEnv)

Check if a simulation has terminated (i.e., has reached the maximum number of periods).
"""
is_terminated(env::SupplyChainEnv) = env.period == env.num_periods

"Set the random seed for a simulation."
Random.seed!(env::SupplyChainEnv) = Random.seed!(env.seed)
