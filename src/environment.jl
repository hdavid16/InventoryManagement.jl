abstract type AbstractEnv end

"""
    Supply Chain Simulation Environment Constructor

# Fields:
- `network::MetaDiGraph`: Supply chain network directed graph.
- `markets::Vector`: Vector of market nodes where demand occurs.
- `producers::Vector`: Vector of producer nodes where material transformation occurs.
- `echelons::Dict`: Dictionary with Vector of nodes downstream of each node in the network (including that node).
- `materials::Vector`: Vector with the names of all materials in the system.
- `inventory_on_hand::DataFrame`: Timeseries with on hand inventories @ each node.
- `inventory_level::DataFrame`: Timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)
- `inventory_pipeline::DataFrame`: Timeseries with pipeline inventories on each arc.
- `inventory_position::DataFrame`: Timeseries with inventory positions @ each node (inventory level + placed replenishments).
- `echelon_stock::DataFrame`: Timeseries with echelon inventory position @ each node.
- `demand::DataFrame`: Timeseries with replenishment orders and market demand placed on each arc.
- `orders::DataFrame`: History of all orders received.
- `open_orders::DataFrame`: Temporary table with outstanding orders.
- `fulfillments::DataFrame`: Fulfilment log for system orders.
- `shipments::DataFrame`: Temporary table with active shipments and time to arrival on each arc.
- `profit::DataFrame`: Timeseries with profit @ each node.
- `metrics::DataFrame`: Service metrics (service level and fill rate) for each supplier and material.
- `reward::Float64`: Final reward in the system (used for RL)
- `period::Int`: Current period in the simulation.
- `num_periods::Int`: Number of periods in the simulation.
- `num_orders::Int`: Number of orders in the system.
- `discount::Float64`: Time discount factor (i.e. interest rate).
- `options::Dict`: Simulation options
  - `backlog::Bool`: Indicator if backlogging is allowed.
  - `reallocate::Bool`: Indicator if unfulfilled requests should be reallocated to alternate suppliers.
  - `evaluate_profit::Bool`: Indicator if the profit should be evaluated at each node.
  - `capacitated_inventory::Bool`: Indicator if inventory limits should be enforced.
  - `guaranteed_service::Bool`: Indicator if simulation should force lost sales after service lead time expires.
  - `adjusted_stock::Bool`: Indicator if the inventory position and echelon stocks should account for orders that have been placed, but are not yet due.
- `seed::Int`: Random seed.
"""
mutable struct SupplyChainEnv <: AbstractEnv
    network::MetaDiGraph
    markets::Vector
    producers::Vector
    echelons::Dict
    materials::Vector
    inventory_on_hand::DataFrame
    inventory_level::DataFrame
    inventory_pipeline::DataFrame
    inventory_position::DataFrame
    echelon_stock::DataFrame
    demand::DataFrame
    orders::DataFrame
    open_orders::DataFrame
    fulfillments::DataFrame
    shipments::DataFrame
    profit::DataFrame
    metrics::DataFrame
    reward::Float64
    period::Int
    num_periods::Int
    num_orders::Int
    discount::Float64
    options::Dict
    seed::Int
end

"""
    SupplyChainEnv(
        network::MetaDiGraph, num_periods::Int;
        backlog::Bool=true, reallocate::Bool=false, 
        guaranteed_service::Bool=false, adjusted_stock::Bool=true,
        capacitated_inventory::Bool=true,
        evaluate_profit::Bool=true, discount::Float64=0., seed::Int=0
    )

Create a `SupplyChainEnv` from a directed graph with metadata (`MetaDiGraph`).
"""
function SupplyChainEnv(
    network::MetaDiGraph, num_periods::Int;
    backlog::Bool=true, reallocate::Bool=false, 
    guaranteed_service::Bool=false, adjusted_stock::Bool=true,
    capacitated_inventory::Bool=true,
    evaluate_profit::Bool=true, discount::Float64=0., seed::Int=0
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
    period, reward, num_orders = 0, 0, 0
    num_periods = num_periods
    if guaranteed_service 
        backlog = true #override backlog if guaranteed_service
    end
    options = Dict(
        :backlog => backlog, 
        :reallocate => reallocate, 
        :evaluate_profit => evaluate_profit,
        :capacitated_inventory => capacitated_inventory,
        :guaranteed_service => guaranteed_service,
        :adjusted_stock => adjusted_stock
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
        num_orders,
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
reset!(env::SupplyChainEnv) = 
    SupplyChainEnv(
        env.network, env.num_periods;
        backlog = env.options[:backlog], reallocate = env.options[:reallocate], 
        guaranteed_service = env.options[:guaranteed_service], adjusted_stock = env.options[:adjusted_stock],
        capacitated_inventory = env.options[:capacitated_inventory],
        evaluate_profit = env.options[:evaluate_profit], discount = env.discount, seed = env.seed
    )

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
    inventory_on_hand = DataFrame(
        period = Int[], 
        node = Int[], 
        material = [], 
        level = Float64[], 
        discarded = []
    )
    #initialize inventory on hand
    for n in nodes, p in mats
        init_inv = get_prop(net, n, :initial_inventory)
        push!(inventory_on_hand, (0, n, p, init_inv[p], 0))
    end

    #pipeline inventory 
    inventory_pipeline = DataFrame(
        period = zeros(Int, length(mats)*length(arcs)),
        arc = repeat(arcs, inner = length(mats)),
        material = repeat(mats, outer = length(arcs)),
        level = zeros(length(mats)*length(arcs))
    )

    #inventory position and level
    inventory_position = select(inventory_on_hand, [:period, :node, :material, :level])
    inventory_level = copy(inventory_position)

    #echenlon inventory position
    echelon_stock = DataFrame(
        period = Int[], 
        node = Int[], 
        material = [], 
        level = Float64[]
    )
    #initialize echelon stocks with initial inventory positions
    for n in nodes, p in mats
        if get_prop(net, n, :inventory_capacity)[p] > 0 #only update if that node can store that material
            ech_pos = sum(filter([:node, :material] => (i1, i2) -> i1 in echelons[n] && i2 == p, inventory_position).level)
            push!(echelon_stock, (0, n, p, ech_pos))
        end
    end

    #replenishment and customer sales orders
    demand = DataFrame(
        period = Int[],
        arc = Tuple[],
        material = Any[],
        quantity = Float64[], 
        fulfilled = Float64[],
        lead = Float64[],
        unfulfilled = Float64[],
        reallocated = Any[]
    )

    #all system orders
    orders = DataFrame(
        id = Int[],
        created = Int[],
        due = Int[],
        arc = Tuple[],
        material = [],
        quantity = []
    )

    #outstanding orders
    open_orders = DataFrame(
        id = Int[],
        arc = Tuple[],
        material = [],
        quantity = [],
        due = []
    )

    #order fulfillments
    fulfillments = DataFrame(
        id = Int[],
        time = [],
        supplier = [],
        amount = []
    )

    #material shipments
    shipments = DataFrame(
        arc = [],
        material = [],
        amount = Float64[],
        lead = Float64[]
    )

    #profit
    profit = DataFrame(
        period = Int[],
        value = Float64[],
        node = Int[]
    )

    #metrics
    metrics = DataFrame()

    return inventory_on_hand, inventory_level, inventory_pipeline, inventory_position, echelon_stock, demand, orders, open_orders, fulfillments, shipments, profit, metrics
end
