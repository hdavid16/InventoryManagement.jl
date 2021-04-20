# InventoryManagement

*Discrete-time simulation environment for Inventory Management in Supply Networks.*

## Overview

*InventoryManagement.jl* allows defining a supply network with three main actors:
- `Producers`: Nodes where inventory transformation takes place (e.g., intermediates or final products are produced).
- `Distributors`: Nodes where inventory is stored and distributed (e.g., distribution centers and retailers).
- `Markets`: Associated with each end `distributor` in the network where end-customers place final product orders.

A `SupplyChainEnv` object is created based on system inputs and network structure, which can be used to simulate stochastic demand at the end distribution centers and inventory replenishment decisions throughout the network. The `SupplyChainEnv` can be used in conjunction with [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl) to train a Reinforcement Learning `agent`.

This package generalizes and extends and the inventory management environment available in [OR-Gym](https://github.com/hubbs5/or-gym).

## Dependencies

*InventoryManagement.jl* relies on the following packages:
- [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl): Define supply network structure and specify node- and edge-specific parameters.
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl): Tabulate results and specify most network parameters.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Define probability distributions for the lead times in between nodes and the market demands at the end distributors.

## Sequence of Events

The following sequence of events occurs in each period of the simulation:
1. 

## Model Inputs

### Node-specific

A `DataFrame` is stored in the metadata of each node using the key `:params` with the following fields:
- `"product"::Vector{String}`: product names
- `"init_inventory"::Vector{Float64}`: initial inventories for each product
- `"holding_cost"::Vector{Float64}`: unit holding cost for each product
- `"production_cost"::Vector{Float64}`: unit production cost for each product
- `"production_capacity"::Vector{Float64}`: maximum production capacity for each product (use `Inf` for uncapacitated production)
- `"market_demand"::Vector{Sampleable}`: probability distributions for the market demands for each product
- `"demand_frequency"::Vector{Float64}`: probability that demand will occur (value between `0.0` and `1.0`)
- `"market_price"::Vector{Float64}`: sales price for each product at the end distributor
- `"demand_penalty"::Vector{Float64}`: unit penalty for unsatisfied demand

### Edge-specific

A `DataFrame` is stored in the metadata of each node using the key `:params` with the following fields:
- `"product"::Vector{String}`: product names
- `"sales_price"::Vector{Float64}`: unit sales price for inventory sent on that edge (from supplier to receiver)
- `"transportation_cost"::Vector{Float64}`: unit transportation cost per period for inventory in-transit

A `Univariate Discrete Distribution` is also defined for the lead time on each edge and stored using the key `:lead_time`.

## `SupplyChainEnv`

A `SupplyChainEnv` has the following fields:
- `network::MetaDiGraph`: Supply Chain Network (metagraph)
- `markets::Array`: list of markets (end distributors)
- `producers::Array`: list of producer nodes
- `distributors::Array`: list of distribution centers (excludes end distributors where markets exist)
- `products::Array`: list of product names (strings)
- `inv_on_hand::DataFrame`: timeseries On Hand Inventory @ each node
- `inv_pipeline::DataFramet`: timeseries Pipeline Inventory on each arc
- `inv_position::DataFrame`: timeseries Inventory Position for each node
- `replenishments::DataFrame`: timeseries Replenishment orders placed on each arc
- `shipments::DataFrame`: current shipments and time to arrival for each node
- `demand::DataFrame`: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market
- `profit::DataFrame`: timeseries with profit at each node
- `reward::Float64`: reward in the system (used for RL)
- `period::Int`: period in the simulation
- `num_periods::Int`: number of periods in the simulation
- `backlog::Bool`: backlogging allowed if true; otherwise, unfulfilled demand is lost sales
- `discount::Float64`: time discount factor (interest rate)

## Example

The example below is for a 30 period simulation of a supply network with two plants (nodes 1 and 2) that supply and end distributor (node 3).

```julia
using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, StatsPlots

#define network connectivity
adjmx = [0 0 1;
         0 0 1;
         0 0 0]
net = MetaDiGraph(adjmx)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_prop!(net, 1, :params,
            DataFrame("product" => ["A", "B", "C", "D", "E"],
                      "init_inventory" => zeros(5),
                      "production_cost" => ones(5),
                      "holding_cost" => zeros(5),
                      "production_capacity" => ones(5)*Inf,
                      "market_demand" => [missing for i in 1:5],
                      "demand_frequency" => [missing for i in 1:5],
                      "demand_penalty" => [missing for i in 1:5],
                      "market_price" => [missing for i in 1:5]))
set_prop!(net, 2, :params,
            DataFrame("product" => ["A", "B", "C", "D", "E"],
                      "init_inventory" => zeros(5),
                      "production_cost" => ones(5),
                      "holding_cost" => ones(5)*0.01,
                      "production_capacity" => ones(5)*Inf,
                      "market_demand" => [missing for i in 1:5],
                      "demand_frequency" => [missing for i in 1:5],
                      "demand_penalty" => [missing for i in 1:5],
                      "market_price" => [missing for i in 1:5]))
set_prop!(net, 3, :params,
            DataFrame("product" => ["A", "B", "C", "D", "E"],
                      "init_inventory" => ones(5)*10,
                      "production_cost" => [missing for i in 1:5],
                      "holding_cost" => ones(5)*0.01,
                      "production_capacity" => [missing for i in 1:5],
                      "market_demand" => [Normal(5,0.5) for i in 1:5],
                      "demand_frequency" => [0.5 for i in 1:5],
                      "demand_penalty" => ones(5)*0.01,
                      "market_price" => ones(5)*3))
#specify sales prices, transportation costs, lead time
set_prop!(net, 1, 3, :params,
            DataFrame("product" => ["A", "B", "C", "D", "E"],
                      "sales_price" => ones(5)*2.00,
                      "transportation_cost" => ones(5)*0.01))
set_prop!(net, 2, 3, :params,
          DataFrame("product" => ["A", "B", "C", "D", "E"],
                    "sales_price" => ones(5)*2.00,
                    "transportation_cost" => ones(5)*0.01))
set_prop!(net, 1, 3, :lead_time, Poisson(5)) #NOTE: could make this specific to MOT
set_prop!(net, 2, 3, :lead_time, Poisson(5)) #NOTE: could make this specific to MOT

# Make environment
env = SupplyChainEnv(net, 30)

# Specify a random replenishment action
action = DataFrame(:product=>["A","B","C","D","E"],
                   [Symbol((a.src,a.dst)) => rand(5)*5 for a in edges(net)]...)

# Run simulation
for t in 1:env.num_periods
    (env)(action)
end

# Plot cumulative profit at each node
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
@df profit plot(:period, :value_cumsum, group=:node)

```
