# InventoryManagement

Discrete-time simulation environment for Inventory Management in Supply Networks.

# Example

```julia
using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, StatsPlots

#define network connectivity
adjmx = [0 0 1;
         0 0 1;
         0 0 0]
net = MetaDiGraph(adjmx)
#define node types (**can be inferred**)
set_prop!(net, 1, :type, :producer)
set_prop!(net, 2, :type, :producer)
set_prop!(net, 3, :type, :distributor)
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

env = SupplyChainEnv(net, 30)
action = DataFrame(:product=>["A","B","C","D","E"],
                   [Symbol((a.src,a.dst)) => rand(5) for a in edges(net)]...)
for t in 1:env.num_periods
    (env)(action)
end

node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
@df profit plot(:period, :value_cumsum, group=:node)

```
