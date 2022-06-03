#2 - echelon system with production
using InventoryManagement

#define network topology and materials:
# Node 1: Plant
#   - Stores raw material B; Converts B => A; receives orders for B from Node 2
#   - Stores product A; receives orders for A from Node 2
# Node 2: Retail; has market demand for A and B
adj_matrix = [1 1;
              0 0]
net = MetaDiGraph(adj_matrix)
set_prop!(net, :materials, [:A, :B])

#specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(
    :initial_inventory => Dict(:B => Inf, :A => 125),
    :holding_cost => Dict(:B => 0.000, :A => 0.002),
    :bill_of_materials => Dict((:B,:A) => -1),
    :supplier_priority => Dict(:A => 1)
))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:A => 100, :B => 50),
    :holding_cost => Dict(:A => 0.003, :B => 0.002),
    :demand_distribution => Dict(:A => 2, :B => 1),
    :demand_frequency => Dict(:A => 1, :B => 1),#Dict(:A => 1/2, :B => 1/3),
    :sales_price => Dict(:A => 3, :B => 2),
    :unfulfilled_penalty => Dict(:A => 0.01, :B => 0.01),
    :supplier_priority => Dict(:A => 1, :B => 1)
))
#specify sales prices, transportation costs, lead time
set_props!(net, 1, 1, Dict(
    :transportation_cost => Dict(:A => 0.1), #production cost
    :lead_time => Dict(:A => 3) #production lead time
))
set_props!(net, 1, 2, Dict(
    :sales_price => Dict(:A => 2, :B => 1),
    :transportation_cost => Dict(:A => 0.01, :B => 0.01),
    :lead_time => Dict(:A => 7, :B => 7))
)

#define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
s = Dict((2,:A) => 50, (2,:B) => 25, (1,:A) => 100) #lower bound on inventory
S = Dict((2,:A) => 100, (2,:B) => 50, (1,:A) => 125) #base stock level

#create environment and run simulation with reorder policy
policy_variable = :inventory_position
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_variable, policy_type, review_period)

#make plots
using DataFrames, StatsPlots
#profit
entity_profit = transform(env.profit, 
    :node => ByRow(i -> i == 1 ? "Plant" : "Retailer") => :entity
)
profit = transform(groupby(entity_profit, :entity), :value => cumsum)
fig1 = @df profit plot(
    :period, :value_cumsum, group=:entity, linetype=:steppost,
    legend = :topleft, xlabel="period", ylabel="cumulative profit"
)

#inventory position
fig2 = @df env.inventory_position plot(
    :period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, 
    legend = :bottomleft, xlabel="period", ylabel="inventory position", yticks = 0:25:125
)
