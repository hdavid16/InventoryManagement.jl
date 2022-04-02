#4-echelon supply chain with market demand at levels 1 and 2 of the chain

using Distributions
using InventoryManagement

##define network connectivity
adj_matrix = [0 1 0 0 0;
              0 0 1 0 0;
              0 0 0 1 0;
              0 0 0 0 1;
              0 0 0 0 0]
net = MetaDiGraph(adj_matrix) # 1 (raw) -> 2 (plant) -> 3 (plant) -> 4 (storage; with demand) -> 5 (retail)
materials = [:A, :B, :C]
bom = [0  0  0; # C -> B -> A
      -1  0  0;
       0 -1  0]
set_prop!(net, :materials, materials)
set_prop!(net, :bill_of_materials, bom)

##specify parameters, holding costs and capacity, market demands and penalty for unfilfilled demand
set_props!(net, 1, Dict(:initial_inventory => Dict(:C => Inf),
                        :inventory_capacity => Dict(:A => 0, :B => 0, :C => Inf)))
set_props!(net, 2, Dict(:initial_inventory => Dict(:C => 215), #initial inventory at plant
                        :inventory_capacity => Dict(:A => 0, :B => 0, :C => Inf),
                        :production_time => Dict(:B => 0)))
set_props!(net, 3, Dict(:initial_inventory => Dict(:B => 190), #initial inventory at plant
                        :inventory_capacity => Dict(:A => 0, :B => Inf, :C => 0),
                        :production_time => Dict(:A => 0)))

set_props!(net, 4, Dict(:initial_inventory => Dict(:A => 105), #initial inventory at storage
                        :inventory_capacity => Dict(:A => Inf, :B => 0, :C => 0),
                        :demand_distribution => Dict(:A => 5), #distribution given
                        :demand_frequency => Dict(:A => 1))) #order every other day on average

set_props!(net, 5, Dict(:initial_inventory => Dict(:A => 60), #initial inventory at retail
                        :inventory_capacity => Dict(:A => Inf, :B => 0, :C => 0),
                        :demand_distribution => Dict(:A => 10), #Normal Distribution double that of storage
                        :demand_frequency => Dict(:A => 1))) #order daily

##specify lead times
set_props!(net, 1, 2, Dict(:lead_time => Dict(:C => 20)))
set_props!(net, 2, 3, Dict(:lead_time => Dict(:B => 15)))
set_props!(net, 3, 4, Dict(:lead_time => Dict(:A => 10)))
set_props!(net, 4, 5, Dict(:lead_time => Dict(:A => 5)))

##define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
policy_variable = :ech_position #variable tracked by policy
s = Dict((2,:C) => 215, (3,:B) => 190, (4,:A) => 165, (5,:A) => 60) #lower bound on inventory
S = Dict((2,:C) => 215, (3,:B) => 190, (4,:A) => 165, (5,:A) => 60) #base stock level

##create environment and run simulation with reorder policy
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_type, review_period, policy_variable)

##make plots
using DataFrames, StatsPlots
#inventory level
df = filter(i -> i.node > 1, env.ech_position)
fig2 = @df df plot(:period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, legend = :bottomright,
                    xlabel="period", ylabel="echelon position")