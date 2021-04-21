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

set_props!(net, 3, Dict(:init_inventory => Dict(p => 10 for p in products),
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

env = SupplyChainEnv(net, 30)

action = DataFrame(:product=>products,
                   [Symbol((a.src,a.dst)) => rand(5)*5 for a in edges(env.network)]...)

for t in 1:env.num_periods
    (env)(action)
end

node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
@df profit plot(:period, :value_cumsum, group=:node)

# action = DataFrame(:product=>products,
#                    [Symbol((a.src,a.dst)) => rand(5)*5 for a in edges(env.network)]...,
#                    [Symbol(n) =>
#                         [rand(Uniform(0,get_prop(env.network, n, :production_capacity)[p]))
#                             for p in products]
#                         for n in env.producers]...)
