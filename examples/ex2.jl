using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, Plots

#define network connectivity
adjmx = [0 0 1;
         0 0 1;
         0 0 0]
net = MetaDiGraph(adjmx)
products = [:A]
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
