using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, StatsPlots

#define network connectivity
net = MetaDiGraph(path_digraph(2)) # 1 -> 2
products = [:A]
set_prop!(net, :products, products)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:production_cost => Dict(p => 0.01 for p in products),
                        :production_time => Dict(p => 0 for p in products),
                        :production_capacity => Dict(p => Inf for p in products)))

set_props!(net, 2, Dict(:initial_inventory => Dict(p => 100 for p in products),
                        :holding_cost => Dict(p => 0.01 for p in products),
                        :demand_distribution => Dict(p => Normal(5,0.5) for p in products),
                        :demand_frequency => Dict(p => 0.5 for p in products),
                        :sales_price => Dict(p => 3 for p in products),
                        :demand_penalty => Dict(p => 0.01 for p in products)))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 2, Dict(:sales_price => Dict(p => 2 for p in products),
                          :transportation_cost => Dict(p => 0.01 for p in products),
                          :lead_time => Poisson(5)))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods)

#define reorder policy parameters
policy = :sS #(s, S) policy
on = :position #monitor inventory position
s = Dict((2,:A) => 20) #lower bound on inventory
S = Dict((2,:A) => 100) #base stock level

#run simulation with reorder policy
for t in 1:env.num_periods
    action = reorder_policy(env, s, S, on, policy)
    (env)(action)
end

#make plots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(:period, :value_cumsum, group=:node, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")

#inventory position
inv_position = filter(i -> i.node in union(env.distributors, env.markets), env.inv_position)
fig2 = @df inv_position plot(:period, :level, group=(:node, :product), linetype=:steppost,
                    xlabel="period", ylabel="inventory position")
