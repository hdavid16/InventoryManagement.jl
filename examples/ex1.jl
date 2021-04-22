using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, StatsPlots

#define network connectivity
net = MetaDiGraph(path_digraph(2)) # 1 -> 2
products = [:A]
set_prop!(net, :products, products)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:production_cost => Dict(:A => 0.01),
                        :production_time => Dict(:A => 0),
                        :production_capacity => Dict(:A => Inf)))

set_props!(net, 2, Dict(:initial_inventory => Dict(:A => 100),
                        :holding_cost => Dict(:A => 0.01),
                        :demand_distribution => Dict(:A => Normal(5,0.5)),
                        :demand_frequency => Dict(:A => 0.5),
                        :sales_price => Dict(:A => 3),
                        :demand_penalty => Dict(:A => 0.01),
                        :supplier_priority => Dict(:A => [1])))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 2, Dict(:sales_price => Dict(:A => 2),
                          :transportation_cost => Dict(:A => 0.01),
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
    action = reorder_policy(env, s, S, on, policy, :priority)
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
