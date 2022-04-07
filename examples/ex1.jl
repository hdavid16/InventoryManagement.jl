#2 - echelon system with production

using Distributions
using InventoryManagement

#define network topology and materials:
# Nodes 1-2: Plant
#   - Node 1: Stores raw material B; Converts B => A; receives orders for B from Node 3
#   - Node 2: Stores product A; receives orders for A from Node 3
#   - Node 3: Retail; has market demand for A and B
adj_matrix = [0 1 1;
              0 0 1;
              0 0 0]
net = MetaDiGraph(adj_matrix)
set_prop!(net, :materials, [:A, :B])

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(
    :initial_inventory => Dict(:B => 150),
    :holding_cost => Dict(:B => 0.01),
    :bill_of_materials => Dict((:B,:A) => -1)
))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:A => 125),
    :holding_cost => Dict(:A => 0.02)
))
set_props!(net, 3, Dict(
    :initial_inventory => Dict(:A => 100, :B => 50),
    :holding_cost => Dict(:A => 0.03, :B => 0.02),
    :demand_distribution => Dict(:A => Normal(3,0.3), :B => Normal(2,0.2)),
    :demand_period => Dict(:A => 2, :B => 3),
    :sales_price => Dict(:A => 3, :B => 2),
    :stockout_penalty => Dict(:A => 0.01, :B => 0.01),
    :supplier_priority => Dict(:A => 2, :B => 1)
))
#specify sales prices, transportation costs, lead time
set_props!(net, 1, 2, Dict(
    :transportation_cost => Dict(:A => 0.1), #production cost
    :lead_time => Dict(:A => 3) #production lead time
))
set_props!(net, 1, 3, Dict(
    :sales_price => Dict(:B => 1),
    :transportation_cost => Dict(:B => 0.1),
    :lead_time => Dict(:B => Poisson(7)))
)
set_props!(net, 2, 3, Dict(
    :sales_price => Dict(:A => 2),
    :transportation_cost => Dict(:A => 0.1),
    :lead_time => Dict(:A => Poisson(7)))
)

#define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
s = Dict((3,:A) => 50, (3,:B) => 25, (2,:A) => 100) #lower bound on inventory
S = Dict((3,:A) => 100, (3,:B) => 50, (2,:A) => 125) #base stock level

#create environment and run simulation with reorder policy
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_type, review_period)

#make plots
using DataFrames, StatsPlots
#profit
node_profit = groupby(env.profit, :node)
profit = transform(node_profit, :value => cumsum)
fig1 = @df profit plot(
    :period, :value_cumsum, group={Node = :node}, 
    legend = :topleft, xlabel="period", ylabel="cumulative profit"
)

#inventory position
inventory_position = filter(i -> i.node == 1 ? i.material == :B : i.node == 2 ? i.material == :A : true, env.inventory_position)
fig2 = @df inventory_position plot(
    :period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, 
    legend = :bottomleft, xlabel="period", ylabel="inventory position", yticks = 0:25:125
)
