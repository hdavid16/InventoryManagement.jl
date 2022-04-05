#2 - echelon system with production

using Distributions
using InventoryManagement

#define network connectivity
adj_matrix = [0 0 1;
              0 0 1;
              0 0 0]
net = MetaDiGraph(adj_matrix) # 1 (plant) & 2 (distribution center) supply 3 (market)
materials = [:A, :B]
bom = [0 0; # B -> A
      -1 0]
set_prop!(net, :materials, materials)
set_prop!(net, :bill_of_materials, bom)

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:initial_inventory => Dict(:A => 75, :B => 100),
                        :holding_cost => Dict(:A => 0, :B => 0),
                        :production_cost => Dict(:A => 0.01)))

set_props!(net, 2, Dict(:initial_inventory => Dict(:A => 125),
                        :holding_cost => Dict(:A => 0)))

set_props!(net, 3, Dict(:initial_inventory => Dict(:A => 100),
                        :holding_cost => Dict(:A => 0.01),
                        :demand_distribution => Dict(:A => Normal(5,0.5)),
                        :demand_frequency => Dict(:A => 2),
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

#define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
s = Dict((3,:A) => 50) #lower bound on inventory
S = Dict((3,:A) => 100) #base stock level

#create environment and run simulation with reorder policy
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_type, review_period)

#make plots
using DataFrames, StatsPlots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(:period, :value_cumsum, group={Node = :node}, legend = :topleft,
                    xlabel="period", ylabel="cumulative profit")

#inventory position
inventory_position = filter(i -> i.node > 1 ? i.material == :A : i.material in [:A,:B], env.inventory_position)
fig2 = @df inventory_position plot(:period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, legend = :bottomleft,
                    xlabel="period", ylabel="inventory position", yticks = 0:25:125)
