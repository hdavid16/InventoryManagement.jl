abstract type AbstractEnv end

"""
    Supply Chain Simulation Environment Constructor

# Fields:
- `network::MetaDiGraph`: Supply chain network directed graph.
- `markets::Vector`: Vector of market nodes where demand occurs.
- `producers::Vector`: Vector of producer nodes where material transformation occurs.
- `echelons::Dict`: Dictionary with Vector of nodes downstream of each node in the network (including that node).
- `materials::Vector`: Vector with the names of all materials in the system.
- `inv_on_hand::DataFrame`: Timeseries with on hand inventories @ each node.
- `inv_level::DataFrame`: Timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)
- `inv_pipeline::DataFrame`: Timeseries with pipeline inventories on each arc.
- `inv_position::DataFrame`: Timeseries with inventory positions @ each node (inventory level + placed replenishments).
- `ech_position::DataFrame`: Timeseries with echelon inventory position @ each node.
- `replenishments::DataFrame`: Timeseries with replenishment orders placed on each arc.
- `shipments::DataFrame`: Temp table with active shipments and time to arrival on each arc.
- `demand::DataFrame`: Timeseries with realization of demand at each market, and amounts sold, unfulfilled demand, and backlog.
- `profit::DataFrame`: Timeseries with profit @ each node.
- `reward::Float64`: Final reward in the system (used for RL)
- `period::Int`: Current period in the simulation.
- `num_periods::Int`: Number of periods in the simulation.
- `discount::Float64`: Time discount factor (i.e. interest rate).
- `options::Dict`: Simulation options
  - `backlog::Bool`: Indicator if backlogging is allowed.
  - `reallocate::Bool`: Indicator if unfulfilled requests should be reallocated to alternate suppliers.
  - `evaluate_profit::Bool`: Indicator if the profit should be evaluated at each node
  - `capacitated_inventory::Bool`: Indicator if inventory limits should be enforced
- `seed::Int`: Random seed.
"""
mutable struct SupplyChainEnv <: AbstractEnv
    network::MetaDiGraph
    markets::Vector
    producers::Vector
    echelons::Dict
    materials::Vector
    inv_on_hand::DataFrame
    inv_level::DataFrame
    inv_pipeline::DataFrame
    inv_position::DataFrame
    ech_position::DataFrame
    replenishments::DataFrame
    shipments::DataFrame
    demand::DataFrame
    profit::DataFrame
    reward::Float64
    period::Int
    num_periods::Int
    discount::Float64
    options::Dict
    seed::Int
end

"""
    SupplyChainEnv(
        network::MetaDiGraph, num_periods::Int;
        discount::Float64=0.0, backlog::Bool=false,
        reallocate::Bool=false, evaluate_profit::Bool=true,
        capacitated_inventory::Bool=true, seed::Int=0
    )

Create a `SupplyChainEnv` from a directed graph with metadata (`MetaDiGraph`).
"""
function SupplyChainEnv(
    network::MetaDiGraph, num_periods::Int;
    discount::Float64=0.0, backlog::Bool=false,
    reallocate::Bool=false, evaluate_profit::Bool=true,
    capacitated_inventory::Bool=true, seed::Int=0
)
    #copy network (avoids issues when changing say num_periods after the Env was already created)
    net = copy(network)
    #get model parameters
    nodes = vertices(net) #network nodes
    echelons = Dict(n => identify_echelons(net, n) for n in nodes) #get nodes in each echelon
    arcs = [(e.src,e.dst) for e in edges(net)] #network edges as Tuples (not Edges)
    mats = get_prop(net, :materials) #materials
    mrkts, plants = identify_nodes(net) #identify nodes 
    #check inputs
    check_inputs!(net, nodes, arcs, mrkts, plants, mats, num_periods)
    #create logging dataframes
    logging_dfs = create_logging_dfs(net, nodes, arcs, mats, echelons)
    #initialize other params
    period, reward = 0, 0
    num_periods = num_periods
    options = Dict(
        :backlog => backlog, 
        :reallocate => reallocate, 
        :evaluate_profit => evaluate_profit,
        :capacitated_inventory => capacitated_inventory
    )
    #create environment
    env = SupplyChainEnv(
        net, 
        mrkts, 
        plants, 
        echelons, 
        mats, 
        logging_dfs...,
        reward, 
        period, 
        num_periods, 
        discount, 
        options, 
        seed
    )
    #seed randomness
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
    filter!(i -> i.period == 0, env.ech_position)
    filter!(i -> i.period == 0, env.replenishments)
    filter!(i -> i.period == 0, env.demand)
    filter!(i -> i.period == 0, env.profit)
    env.shipments = DataFrame(
        arc = [],
        material = [],
        amount = Float64[],
        lead = Int[]
    )

    Random.seed!(env)
end

"""
    is_terminated(env::SupplyChainEnv)

Check if a simulation has terminated (i.e., has reached the maximum number of periods).
"""
is_terminated(env::SupplyChainEnv) = env.period == env.num_periods

"Set the random seed for a simulation."
Random.seed!(env::SupplyChainEnv) = Random.seed!(env.seed)

"""
    create_logging_dfs(net::MetaDiGraph, nodes::Vector, arcs::Vector, mats::Vector, echelons::Dict)

Create `DataFrames` to log timeseries results.
"""
function create_logging_dfs(net::MetaDiGraph, nodes::Base.OneTo, arcs::Vector, mats::Vector, echelons::Dict)
    #inventory on hand
    inv_on_hand = DataFrame(
        period = Int[], 
        node = Int[], 
        material = [], 
        level = Float64[], 
        discarded = []
    )
    #initialize inventory on hand
    for n in nodes, p in mats
        init_inv = get_prop(net, n, :initial_inventory)
        push!(inv_on_hand, (0, n, p, init_inv[p], 0))
    end

    #pipeline inventory 
    inv_pipeline = DataFrame(
        period = zeros(Int, length(mats)*length(arcs)),
        arc = repeat(arcs, inner = length(mats)),
        material = repeat(mats, outer = length(arcs)),
        level = zeros(length(mats)*length(arcs))
    )

    #inventory position and level
    inv_position = select(inv_on_hand, [:period, :node, :material, :level])
    inv_level = copy(inv_position)

    #echenlon inventory position
    ech_position = DataFrame(
        period = Int[], 
        node = Int[], 
        material = [], 
        level = Float64[]
    )
    #initialize echelon positions with initial inventory positions
    for n in nodes, p in mats
        if get_prop(net, n, :inventory_capacity)[p] > 0 #only update if that node can store that material
            ech_pos = sum(filter([:node, :material] => (i1, i2) -> i1 in echelons[n] && i2 == p, inv_position).level)
            push!(ech_position, (0, n, p, ech_pos))
        end
    end

    #replenishment orders
    replenishments = DataFrame(
        period = Int[],
        arc = Tuple[],
        material = Any[],
        requested = Float64[],
        accepted = Float64[],
        lead = Float64[],
        unfulfilled = Float64[],
        reallocated = Any[]
    )

    #material shipments
    shipments = DataFrame(
        arc = [],
        material = [],
        amount = Float64[],
        lead = Float64[]
    )

    #finished good demands
    demand = DataFrame(
        period = Int[],
        node = Int[],
        material = Any[],
        demand = Float64[],
        sold = Float64[],
        unfulfilled = Float64[]
    )

    #profit
    profit = DataFrame(
        period = Int[],
        value = Float64[],
        node = Int[]
    )

    return inv_on_hand, inv_level, inv_pipeline, inv_position, ech_position, replenishments, shipments, demand, profit
end
