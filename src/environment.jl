abstract type AbstractEnv end

mutable struct SupplyChainEnv <: AbstractEnv
    network::MetaDiGraph #Supply Chain Network
    markets::Array #Array of markets (end distributors)
    producers::Array #Array of producer nodes
    distributors::Array #Array of distribution centers (excludes market nodes)
    products::Array #Array of product names (strings)
    inv_on_hand::DataFrame #Timeseries On Hand Inventory @ each node
    inv_pipeline::DataFrame #Timeseries Pipeline Inventory on each arc
    inv_position::DataFrame #Timeseries Inventory Position for each node
    replenishments::DataFrame #Timeseries Replenishment orders placed on each arc
    shipments::DataFrame #Current shipments and time to arrival for each node
    production::DataFrame #Current material production commited to an arc and time to ship
    demand::DataFrame #Timeseries with realization of demand at each market
    profit::DataFrame #Timeseries with profit at each node
    reward::Float64 #Final reward in the system (used for RL)
    period::Int #period in the simulation
    num_periods::Int #number of periods in the simulation
    backlog::Bool #Backlogging allowed if true; otherwise, lost sales
    discount::Float64 #Time discount factor
    seed::Int #Random seed
end

function SupplyChainEnv(network::MetaDiGraph, num_periods::Int;
                        backlog::Bool=true, discount::Float64=0.0,
                        seed::Int=0)
    #perform checks
        #order of products in init_inventory, production_capacity, market_demand and demand_frequency keys must be the same

    #get main nodes
    nodes = vertices(network)
    #get main edges
    arcs = [(e.src,e.dst) for e in edges(network)]
    #get end distributors, producers, and distribution centers
    mrkts = [n for n in nodes if isempty(outneighbors(network,n))]
    plants = [n for n in nodes if isempty(inneighbors(network,n))] #assumes plants are at the top
    dcs = setdiff(nodes, mrkts, plants)
    #get products
    prods = get_prop(network, :products)
    #create logging dataframes
    inv_on_hand = DataFrame(:period => Int[], :node => Int[], :product => [], :level => Float64[])
    for n in setdiff(nodes, plants), p in prods
        init_inv = get_prop(network, n, :init_inventory)
        push!(inv_on_hand, (0, n, p, init_inv[p]))
    end
    inv_pipeline = DataFrame(:period => zeros(Int, length(prods)*length(arcs)),
                             :arc => repeat(arcs, inner = length(prods)),
                             :product => repeat(prods, outer = length(arcs)),
                             :level => zeros(length(prods)*length(arcs)))
    inv_position = copy(inv_on_hand)
    replenishments = DataFrame(:period => zeros(Int, length(prods)*length(arcs)),
                               :arc => repeat(arcs, inner = length(prods)),
                               :product => repeat(prods, outer = length(arcs)),
                               :amount => zeros(length(prods)*length(arcs)),
                               :lead => zeros(Int, length(prods)*length(arcs)))
    shipments = DataFrame(:arc => [],
                          :product => [],
                          :amount => Float64[],
                          :lead => Int[])
    production = DataFrame(:arc => [],
                           :product => [],
                           :amount => Float64[],
                           :lead => Int[])
    demand = DataFrame(:period => zeros(Int, length(prods)*length(mrkts)),
                       :node => repeat(mrkts, inner = length(prods)),
                       :product => repeat(prods, outer = length(mrkts)),
                       :demand => zeros(length(prods)*length(mrkts)),
                       :sale => zeros(length(prods)*length(mrkts)),
                       :unfulfilled => zeros(length(prods)*length(mrkts)),
                       :backlog => zeros(length(prods)*length(mrkts)))
    profit = DataFrame(:period => zeros(Int, length(nodes)),
                       :value => zeros(length(nodes)),
                       :node => nodes)
    reward = 0
    period = 0
    num_periods = num_periods
    SupplyChainEnv(network, mrkts, plants, dcs, prods, inv_on_hand, inv_pipeline, inv_position,
                    replenishments, shipments, production, demand,
                    profit, reward, period, num_periods, backlog, discount, seed)
end

Random.seed!(env::SupplyChainEnv) = Random.seed!(env.seed)
