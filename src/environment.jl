abstract type AbstractEnv end

"""
    Supply Chain Simulation Environment Constructor

# Fields:
- `network::MetaDiGraph`: Supply chain network directed graph.
- `markets::Vector`: Vector of market nodes where demand occurs.
- `producers::Vector`: Vector of producer nodes where material transformation occurs.
- `echelons::Dict`: Dictionary with Vector of nodes downstream of each node in the network (including that node).
- `materials::Vector`: Vector with the names of all materials in the system.
- `tmp::Dict`: Dictionary storing the current and previous on-hand, pipeline, echelon, discarded inventories for each material at each node.
- `inventory::DataFrame`: Timeseries with inventories @ each node.
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
    tmp::Dict
    inventory::DataFrame
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
    mats = get_prop(net, :materials) #materials
    mrkts, plants = identify_nodes(net) #identify nodes 
    #check inputs
    check_inputs!(net, mrkts, plants)
    #create current and previous inventory placeholders
    tmp = create_temp_placeholders(net, echelons)
    #create logging dataframes
    logging_dfs = create_logging_dfs(net, tmp)
    #initialize other params
    period, reward, num_orders = 0, 0, 0
    num_periods = num_periods
    options = Dict(
        :backlog => guaranteed_service ? true : backlog, #override backlog if guaranteed_service
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
        tmp,
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
    create_logging_dfs(net::MetaDiGraph, tmp::Dict)

Create `DataFrames` to log timeseries results.
"""
function create_logging_dfs(net::MetaDiGraph, tmp::Dict)
    nodes = vertices(net) #network nodes
    arcs = [(src(e),dst(e)) for e in edges(net)] #network edges as Tuples (not Edges)
    mats = get_prop(net, :materials) #materials

    #inventory
    inventory = DataFrame(
        period = Int[],
        location = Any[],
        material = Any[],
        amount = Real[],
        type = Symbol[]
    )

    #initialize inventories (on-hand, position, level, echelon)
    for n in nodes
        for p in get_prop(net,n,:node_materials)
            push!(inventory, (0, n, p, tmp[n,p,:on_hand], :on_hand))
            push!(inventory, (0, n, p, tmp[n,p,:on_hand], :position))
            push!(inventory, (0, n, p, tmp[n,p,:on_hand], :level))
            push!(inventory, (0, n, p, tmp[n,p,:echelon], :echelon))
        end
    end

    #pipeline inventory 
    for a in arcs
        for p in get_prop(net,a[1],:node_materials)
            push!(inventory, (0, a, p, 0, :pipeline))
        end
    end

    #all system orders
    orders = DataFrame(
        id = Int[],
        created = Int[],
        due = Int[],
        arc = Tuple[],
        material = Any[],
        amount = Real[]
    )

    #outstanding orders
    open_orders = DataFrame(
        id = Int[],
        due = Int[],
        arc = Tuple[],
        material = Any[],
        amount = Any[]
    )

    #order fulfillments
    fulfillments = DataFrame(
        id = Int[],
        sent = Int[],
        delivered = Int[],
        arc = Tuple[],
        material = Any[],
        amount = Any[]
    )

    #material shipments
    shipments = DataFrame(
        arc = Tuple[],
        material = Any[],
        amount = Real[],
        lead = Int[]
    )

    #profit
    profit = DataFrame(
        period = Int[],
        value = Real[],
        node = Int[]
    )

    #metrics
    metrics = DataFrame()

    return inventory, orders, open_orders, fulfillments, shipments, profit, metrics
end

"""
    create_temp_placeholders(net::MetaDiGraph, echelons::Dict)

Create temporary placeholders to store current and previous on-hand, pipeline, echelon, 
    and discarded inventory. This will avoid filtering many times during the simulation.
"""
function create_temp_placeholders(net::MetaDiGraph, echelons::Dict)
    tmp = Dict()
    i0 = Dict(n => get_prop(net, n, :initial_inventory) for n in vertices(net))
    for n in vertices(net)
        for m in get_prop(net, n, :node_materials)
            ech_stk = sum([m in keys(i0[n]) ? i0[n][m] : 0 for n1 in echelons[n]]) #sum downstream inventories for that material
            tmp[n,m,:echelon] = ech_stk
            tmp[n,m,:on_hand] = i0[n][m]
            tmp[n,m,:position] = i0[n][m]
            tmp[n,m,:level] = i0[n][m]
            tmp[n,m,:discarded] = 0
            for dst in outneighbors(net,n)
                tmp[(n,dst),m,:pipeline] = 0
            end
        end
    end

    return tmp
end