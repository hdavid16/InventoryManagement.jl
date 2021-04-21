# InventoryManagement.jl

Discrete-time simulation environment for Inventory Management in Supply Networks.

## Overview

*InventoryManagement.jl* allows defining a supply network with three main actors:
- `Producers`: Nodes where inventory transformation takes place (e.g., intermediates or final products are produced). These are the top-most (source) nodes in the network.
- `Distributors`: Intermediate nodes where inventory is stored and distributed (e.g., distribution centers).
- `Markets`: Nodes where end-customers place final product orders. These are the last  (sink) nodes in the network.

A `SupplyChainEnv` object is created based on system inputs and network structure, which can be used to simulate stochastic demand at the end distribution centers and inventory replenishment decisions throughout the network. The `SupplyChainEnv` can be used in conjunction with [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl) to train a Reinforcement Learning `agent`.

This package generalizes and extends and the inventory management environment available in [OR-Gym](https://github.com/hubbs5/or-gym).

## Dependencies

*InventoryManagement.jl* relies on the following packages:
- [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl): Define supply network structure and specify node- and edge-specific parameters.
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl): Tabulate results and specify most network parameters.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Define probability distributions for the lead times in between nodes and the market demands at the end distributors.

## Sequence of Events

The following sequence of events occurs in each period of the simulation:
1. Start period.
2. Place inventory replenishment orders at each node. These are limited to the available production capacity (supplier is a producer) or available inventory (supplier is a distribution node).
   - Distributors ship inventory.
   - Producers manufacture products (a production lead time begins). Production costs are incurred at the start of production.
   - Producers send orders that have completed (after the production lead time).
4. Receive inventory that has arrived at each node (after the lead time has transpired).
5. Pay suppliers for inventory received.
6. Pay shipper for inventory shipped.
7. Market demand occurs after tossing a weighted coin with the probability of demand occurring defined by the `demand_frequency`.
8. Demand (including any backlog if `backlog = true`) is fulfilled up to available inventory at the markets.
9. Unfulfilled demand is penalized and backlogged (if `backlog = true`).
10. Each node pays a holding cost and a transportation cost for on-hand inventory and in-transit inventory at each period.

## Model Inputs

### Node-specific

`Producers` will have the following fields in their node metadata:
- `:production_cost::Dict`: unit production cost for each product (`keys`)
- `:production_capacity::Dict`: maximum production capacity for each product (`keys`).
- `:production_time::Dict`: production lead time for each product (`keys`).

`Distributors` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each product (`keys`)
- `:holding_cost::Dict`: unit holding cost for each product (`keys`)

`Markets` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each product (`keys`)
- `:holding_cost::Dict`: unit holding cost for each product (`keys`)
- `:demand_distribution::Dict`: probability distributions for the market demands for each product (`keys`)
- `:demand_frequency::Dict`: probability that demand will occur (value between `0.0` and `1.0`) for each product (`keys`)
- `:sales_price::Dict`: market sales price for each product (`keys`)
- `:demand_penalty::Dict`: unit penalty for unsatisfied market demand for each product (`keys`)

### Edge-specific

All edges have the following fields in their metadata:
- `:sales_price::Dict`: unit sales price for inventory sent on that edge (from supplier to receiver) for each product (`keys`)
- `:transportation_cost::Dict`: unit transportation cost per period for inventory in-transit for each product (`keys`)
- `:lead_time::Distribution{Univariate, Discrete}`: the lead time on each edge

## Model Output

A `SupplyChainEnv` has the following fields:
- `network::MetaDiGraph`: Supply Chain Network (metagraph)
- `markets::Array`: list of market nodes
- `producers::Array`: list of producer nodes
- `distributors::Array`: list of distribution nodes (excludes end distributors where markets exist)
- `products::Array`: list of product names (strings)
- `inv_on_hand::DataFrame`: timeseries On Hand Inventory @ each node at the end of each period
- `inv_pipeline::DataFramet`: timeseries Pipeline Inventory on each edge at the end of each period
- `inv_position::DataFrame`: timeseries Inventory Position for each node at the end of each period
- `replenishments::DataFrame`: timeseries Replenishment orders placed on each edge at the end of each period
- `shipments::DataFrame`: current shipments and time to arrival for each node
- `production::DataFrame`: current material production committed to an edge and lead time to ship
- `demand::DataFrame`: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market
- `profit::DataFrame`: timeseries with profit at each node
- `reward::Float64`: reward in the system (used for RL)
- `period::Int`: period in the simulation
- `num_periods::Int`: number of periods in the simulation
- `backlog::Bool`: backlogging allowed if true; otherwise, unfulfilled demand is lost sales
- `discount::Float64`: time discount factor (interest rate)
- `seed::Int`: random seed

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
products = [:A, :B, :C, :D, :E]
set_prop!(net, :products, products)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:production_cost => Dict(p => 0.01 for p in products),
                        :production_time => Dict(p => 1 for p in products),
                        :production_capacity => Dict(p => Inf for p in products)))

set_props!(net, 2, Dict(:production_cost => Dict(p => 0.01 for p in products),
                        :production_time => Dict(p => 1 for p in products),
                        :production_capacity => Dict(p => Inf for p in products)))

set_props!(net, 3, Dict(:initial_inventory => Dict(p => 10 for p in products),
                        :holding_cost => Dict(p => 0.01 for p in products),
                        :demand_distribution => Dict(p => Normal(5,0.5) for p in products),
                        :demand_frequency => Dict(p => 0.5 for p in products),
                        :sales_price => Dict(p => 3 for p in products),
                        :demand_penalty => Dict(p => 0.01 for p in products)))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 3, Dict(:sales_price => Dict(p => 2 for p in products),
                          :transportation_cost => Dict(p => 0.01 for p in products),
                          :lead_time => Poisson(5)))

set_props!(net, 2, 3, Dict(:sales_price => Dict(p => 2 for p in products),
                          :transportation_cost => Dict(p => 0.01 for p in products),
                          :lead_time => Poisson(5)))

#create environment
env = SupplyChainEnv(net, 30)

#define action
action = DataFrame(:product=>products,
                   [Symbol((a.src,a.dst)) => rand(5)*5 for a in edges(env.network)]...)

#run simulation
for t in 1:env.num_periods
    (env)(action)
end

#make plots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig = @df profit plot(:period, :value_cumsum, group=:node, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")
display(fig)
#on hand inventory
onhand = filter(i -> i.node in union(env.distributors, env.markets), env.inv_on_hand)
fig = @df onhand plot(:period, :level, group=(:node, :product), linetype=:steppost,
                    xlabel="period", ylabel="on hand inventory level")
display(fig)
#inventory position
position = filter(i -> i.node in union(env.distributors, env.markets), env.inv_position)
fig = @df position plot(:period, :level, group=(:node, :product), linetype=:steppost,
                    xlabel="period", ylabel="inventory position")
display(fig)
#production
production = filter(i -> i.arc[1] in env.producers, env.replenishments)
transform!(production, :arc .=> ByRow(y -> y[1]) .=> :plant)
fig = @df production plot(:period, :amount, group=(:plant,:product), linetype=:steppost,
                xlabel="period", ylabel="units produced")
display(fig)
#demand profile
demand = filter(i -> i.node in env.markets, env.demand)
fig = @df demand plot(:period, :demand, group=:product, linetype=:steppost,
                xlabel="period", ylabel="demand")
display(fig)
fig = @df demand plot(:period, :unfulfilled, group=:product, legend=:topleft, linetype=:steppost,
                xlabel="period", ylabel="unfulfilled demand")
display(fig)

```
