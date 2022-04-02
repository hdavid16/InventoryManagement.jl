abstract type AbstractEnv end

"""
    Supply Chain Simulation Environment Constructor

# Fields:
- `network::MetaDiGraph`: Supply chain network directed graph.
- `markets::Vector`: Vector of market nodes where demand occurs.
- `producers::Vector`: Vector of producer nodes where material transformation occurs.
- `distributors::Vector`: Vector of distribution centers (excludes market nodes).
- `echelons::Dict`: Dictionary with Vector of nodes downstream of each node in the network (including that node).
- `materials::Vector`: Vector with the names of all materials in the system.
- `products::Vector`: Vector with the names of all final products in the system.
- `bill_of_materials::Array`: Square matrix with BOM (rows = inputs, cols = ouputs); indices follow materials list; row with positive value is a co-product, row with negative value is a input (reactant).
- `inv_on_hand::DataFrame`: Timeseries with on hand inventories @ each node.
- `inv_level::DataFrame`: Timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)
- `inv_pipeline::DataFrame`: Timeseries with pipeline inventories on each arc.
- `inv_position::DataFrame`: Timeseries with inventory positions @ each node (inventory level + placed replenishments).
- `ech_position::DataFrame`: Timeseries with echelon inventory position @ each node.
- `replenishments::DataFrame`: Timeseries with replenishment orders placed on each arc.
- `shipments::DataFrame`: Temp table with active shipments and time to arrival on each arc.
- `production::DataFrame`: Temp table with active material production commited to an arc and time to ship.
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
    distributors::Vector
    echelons::Dict
    materials::Vector
    products::Vector
    bill_of_materials::Array
    inv_on_hand::DataFrame
    inv_level::DataFrame
    inv_pipeline::DataFrame
    inv_position::DataFrame
    ech_position::DataFrame
    replenishments::DataFrame
    shipments::DataFrame
    production::DataFrame
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
    SupplyChainEnv(network::MetaDiGraph, num_periods::Int;
                        discount::Float64=0.0, backlog::Bool=false,
                        reallocate::Bool=false, evaluate_profit::Bool=true,
                        capacitated_inventory::Bool=true, seed::Int=0)

Create a `SupplyChainEnv` from a directed graph with metadata (`MetaDiGraph`).
"""
function SupplyChainEnv(network::MetaDiGraph, num_periods::Int;
                        discount::Float64=0.0, backlog::Bool=false,
                        reallocate::Bool=false, evaluate_profit::Bool=true,
                        capacitated_inventory::Bool=true, seed::Int=0)
    #copy network (avoids issues when changing say num_periods after the Env was already created)
    net = copy(network)
    #get main nodes
    nodes = vertices(net)
    #get nodes in each echelon
    echelons = Dict(n => identify_echelons(net, n) for n in nodes)
    #get main edges
    arcs = [(e.src,e.dst) for e in edges(net)]
    #get end distributors, producers, and distribution centers
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :demand_penalty, :demand_sequence] #keys to identify a market
    plant_keys = [:production_cost, :production_time, :production_capacity] #keys to identify a plant (producer)
    mrkts = [n for n in nodes if !isempty(intersect(market_keys, keys(net.vprops[n])))]
    plants = [n for n in nodes if !isempty(intersect(plant_keys, keys(net.vprops[n])))]
    sinks = [n for n in nodes if isempty(outneighbors(net, n))]
    dcs = setdiff(nodes, plants, sinks)
    #get materials
    mats = get_prop(net, :materials)
    prods = union(
        [
            intersect([:demand_distribution, :demand_sequence], keys(net.vprops[n]))[1] |> #find which of the keys is used in the market node (demand_distribution or demand_sequence)
            param -> findall(i -> i isa Sampleable || !iszero(i), get_prop(net,n,param)) 
            for n in mrkts
        ]...
    )
    #check inputs
    check_inputs!(net, nodes, arcs, mrkts, plants, mats, num_periods)
    #get bill of materials
    bom = get_prop(net, :bill_of_materials)
    #get material conversion to products
    _, conversion_dict = material_conversion(net)
    set_prop!(net, :conversion_dictionary, conversion_dict)
    #create logging dataframes
    inv_on_hand = DataFrame(:period => Int[], :node => Int[], :material => [], :level => Float64[], :discarded => [])
    for n in nodes, p in mats
        init_inv = get_prop(net, n, :initial_inventory)
        push!(inv_on_hand, (0, n, p, init_inv[p], 0))
    end
    inv_pipeline = DataFrame(:period => zeros(Int, length(mats)*length(arcs)),
                             :arc => repeat(arcs, inner = length(mats)),
                             :material => repeat(mats, outer = length(arcs)),
                             :level => zeros(length(mats)*length(arcs)))
    inv_position = select(inv_on_hand, [:period, :node, :material, :level])
    inv_level = copy(inv_position)
    ech_position = DataFrame(:period => Int[], :node => Int[], :material => [], :level => Float64[])
    for n in nodes, p in mats
        if get_prop(net, n, :inventory_capacity)[p] > 0 #only update if that node can store that material
            ech_pos = sum(filter([:node, :material] => (i1, i2) -> i1 in echelons[n] && i2 == p, inv_position).level)
            push!(ech_position, (0, n, p, ech_pos))
        end
    end
    replenishments = DataFrame(:period => Int[],#zeros(Int, length(mats)*length(arcs)),
                               :arc => Tuple[],#repeat(arcs, inner = length(mats)),
                               :material => Any[],#repeat(mats, outer = length(arcs)),
                               :requested => Float64[],#total amount requested
                               :accepted => Float64[],#zeros(length(mats)*length(arcs)),
                               :lead => Float64[],#zeros(Int, length(mats)*length(arcs)),
                               :unfulfilled => Float64[],#zeros(length(mats)*length(arcs)),
                               :reallocated => Any[])#Any[missing for i in 1:length(mats)*length(arcs)])
    shipments = DataFrame(:arc => [],
                          :material => [],
                          :amount => Float64[],
                          :lead => Float64[])
    production = DataFrame(:arc => [],
                           :material => [],
                           :amount => Float64[],
                           :lead => Float64[])
    demand = DataFrame(:period => Int[],#zeros(Int, length(mats)*length(mrkts)),
                       :node => Int[],#repeat(mrkts, inner = length(mats)),
                       :material => Any[],#repeat(mats, outer = length(mrkts)),
                       :demand => Float64[],#zeros(length(mats)*length(mrkts)),
                       :sold => Float64[],#zeros(length(mats)*length(mrkts)),
                       :unfulfilled => Float64[])#zeros(length(mats)*length(mrkts)))
    profit = DataFrame(:period => Int[],#zeros(Int, length(nodes)),
                       :value => Float64[],#zeros(length(nodes)),
                       :node => Int[])#nodes)
    reward = 0
    period = 0
    num_periods = num_periods
    options = Dict(:backlog => backlog, :reallocate => reallocate, 
                   :evaluate_profit => evaluate_profit,
                   :capacitated_inventory => capacitated_inventory)
    env = SupplyChainEnv(
        net, mrkts, plants, dcs, echelons, mats, prods, bom, 
        inv_on_hand, inv_level, inv_pipeline, inv_position, ech_position, 
        replenishments, shipments, production, demand,
        profit, reward, period, num_periods, discount, options, seed
    )

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
