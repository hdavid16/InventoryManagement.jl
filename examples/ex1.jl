using LightGraphs, MetaGraphs, DataFrames, Distributions
using InventoryManagement, StatsPlots

#define network connectivity
net = MetaDiGraph(path_digraph(2)) # 1 -> 2
products = [:A, :B]
bom = [0 0; # B -> A
      -1 0]
set_prop!(net, :products, products)
set_prop!(net, :bill_of_materials, bom)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:initial_inventory => Dict(:A => 0, :B => 100),
                        :holding_cost => Dict(:A => 0, :B => 0),
                        :production_cost => Dict(:A => 0.01, :B => 0),
                        :production_time => Dict(:A => 0, :B => 0),
                        :production_capacity => Dict(:A => Inf, :B => 0)))

set_props!(net, 2, Dict(:initial_inventory => Dict(:A => 100, :B => 0),
                        :holding_cost => Dict(:A => 0.01, :B => 0),
                        :demand_distribution => Dict(:A => Normal(5,0.5),
                                                     :B => zeros(2)),
                        :demand_frequency => Dict(:A => 0.5, :B => 0),
                        :sales_price => Dict(:A => 3, :B => 0),
                        :demand_penalty => Dict(:A => 0.01, :B => 0),
                        :supplier_priority => Dict(:A => [1], :B => [1])))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 2, Dict(:sales_price => Dict(:A => 2, :B => 0),
                          :transportation_cost => Dict(:A => 0.01, :B => 0),
                          :lead_time => Poisson(5)))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods)

#define reorder policy parameters
policy = :sS #(s, S) policy
on = :position #monitor inventory position
s = Dict((2,:A) => 20, (2,:B) => 0) #lower bound on inventory
S = Dict((2,:A) => 100, (2,:B) => 0) #base stock level

#run simulation with reorder policy
for t in 1:env.num_periods
    action = reorder_policy(env, s, S, on, policy, :priority)
    (env)(action)
end

#make plots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(:period, :value_cumsum, group={Node = :node}, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")

#inventory position
fig2 = @df env.inv_position plot(:period, :level, group={Node = :node, Product = :product}, linetype=:steppost,
                    xlabel="period", ylabel="inventory position")
