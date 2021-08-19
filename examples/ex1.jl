using LightGraphs, MetaGraphs, Distributions
using InventoryManagement

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
                        :production_cost => Dict(:A => 0.01),
                        :production_time => Dict(:A => 0),
                        :production_capacity => Dict(:A => Inf)))

set_props!(net, 2, Dict(:initial_inventory => Dict(:A => 125),
                        :inventory_capacity => Dict(:A => Inf),
                        :holding_cost => Dict(:A => 0)))

set_props!(net, 3, Dict(:initial_inventory => Dict(:A => 100),
                        :inventory_capacity => Dict(:A => Inf),
                        :holding_cost => Dict(:A => 0.01),
                        :demand_distribution => Dict(:A => Normal(5,0.5)),
                        :demand_frequency => Dict(:A => 0.5),
                        :sales_price => Dict(:A => 3),
                        :demand_penalty => Dict(:A => 0.01),
                        :supplier_priority => Dict(:A => [1,2])))

#specify sales prices, transportation costs, lead time
set_props!(net, 1, 3, Dict(:sales_price => Dict(:A => 2),
                          :transportation_cost => Dict(:A => 0.1),
                          :lead_time => Dict(:A => Poisson(3))))

set_props!(net, 2, 3, Dict(:sales_price => Dict(:A => 1),
                          :transportation_cost => Dict(:A => 0.1),
                          :lead_time => Dict(:A => Poisson(7))))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)

#define reorder policy parameters
policy = :sS #(s, S) policy
freq = 1 #continuous review
s = Dict((3,:A) => 50) #lower bound on inventory
S = Dict((3,:A) => 100) #base stock level

#run simulation with reorder policy
simulate_policy!(env, s, S, policy, freq)

#make plots
using DataFrames, StatsPlots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(:period, :value_cumsum, group={Node = :node}, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")

#inventory position
inv_position = filter(i -> i.node > 1 ? i.material == :A : i.material in [:A,:B], env.inv_position)
fig2 = @df inv_position plot(:period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, legend = :bottomleft,
                    xlabel="period", ylabel="inventory position", yticks = 0:25:125)
