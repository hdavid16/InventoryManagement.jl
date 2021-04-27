using LightGraphs, MetaGraphs, Distributions
using InventoryManagement, StatsPlots, Random

#define network connectivity
adj_matrix = [0 0 1;
              0 0 1;
              0 0 0]
net = MetaDiGraph(adj_matrix) # 1 & 2 supply 3
materials = [:A, :B]
bom = [0 0; # B -> A
      -1 0]
set_prop!(net, :materials, materials)
set_prop!(net, :bill_of_materials, bom)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:initial_inventory => Dict(:A => 75, :B => 100),
                        :inventory_capacity => Dict(:A => Inf, :B => Inf),
                        :holding_cost => Dict(:A => 0, :B => 0),
                        :production_cost => Dict(:A => 0.01, :B => 0),
                        :production_time => Dict(:A => 0, :B => 0),
                        :production_capacity => Dict(:A => Inf, :B => 0)))

set_props!(net, 2, Dict(:initial_inventory => Dict(:A => 125, :B => 0),
                        :inventory_capacity => Dict(:A => Inf, :B => Inf),
                        :holding_cost => Dict(:A => 0, :B => 0)))

set_props!(net, 3, Dict(:initial_inventory => Dict(:A => 100, :B => 0),
                        :inventory_capacity => Dict(:A => Inf, :B => Inf),
                        :holding_cost => Dict(:A => 0.01, :B => 0),
                        :demand_distribution => Dict(:A => Normal(5,0.5),
                                                     :B => [0]),
                        :demand_frequency => Dict(:A => 0.5, :B => 0),
                        :sales_price => Dict(:A => 3, :B => 0),
                        :demand_penalty => Dict(:A => 0.01, :B => 0),
                        :supplier_priority => Dict(:A => [1,2], :B => [1,2])))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 3, Dict(:sales_price => Dict(:A => 2, :B => 0),
                          :transportation_cost => Dict(:A => 0.01, :B => 0),
                          :lead_time => Poisson(3)))

set_props!(net, 2, 3, Dict(:sales_price => Dict(:A => 1, :B => 0),
                          :transportation_cost => Dict(:A => 0.01, :B => 0),
                          :lead_time => Poisson(7)))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods)
Random.seed!(env) #set random seed

#define reorder policy parameters
policy = :sS #(s, S) policy
on = :position #monitor inventory position
s = Dict((3,:A) => 50, (3,:B) => 0) #lower bound on inventory
S = Dict((3,:A) => 100, (3,:B) => 0) #base stock level

#run simulation with reorder policy
simulate_policy!(env, s, S, on, policy, :priority)

#make plots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(:period, :value_cumsum, group={Node = :node}, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")

#inventory position
inv_position = filter(i -> i.node > 1 ? i.material == :A : i.material in [:A,:B], env.inv_position)
fig2 = @df inv_position plot(:period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, legend = :bottomright,
                    xlabel="period", ylabel="inventory position", yticks = 0:25:125)
