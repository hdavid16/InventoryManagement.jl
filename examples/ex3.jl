#3-echelon supply chain with market demand at levels 1 and 2 of the chain
using Distributions
using InventoryManagement

##define network topology and materials involved:
# Node 1 = Raw material supplier
# Nodes 2-4 = Plant
#     - Node 2 = Store raw C; Convert C => B
#     - Node 3 = Store intermediate B; Convert B => A
#     - Node 4 = Store product A (direct demand of A occurs here)
# Node 5 = Retail (direct demand of A occurs here)
net = MetaDiGraph(5)
connect_nodes!(net,
    1 => 2,
    2 => 3,
    3 => 4,
    4 => 5
)
set_prop!(net, :materials, [:A, :B, :C])

##specify node parameters
set_props!(net, 1, Dict(
    :initial_inventory => Dict(:C => Inf),
    :inventory_capacity => Dict(:A => 0, :B => 0, :C => Inf)
))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:C => 215), #initial inventory at plant
    :inventory_capacity => Dict(:A => 0, :B => 0, :C => Inf),
    :bill_of_materials => Dict((:C,:B) => -1) #C => B
))
set_props!(net, 3, Dict(
    :initial_inventory => Dict(:B => 190), #initial inventory at plant
    :inventory_capacity => Dict(:A => 0, :B => Inf, :C => 0),
    :bill_of_materials => Dict((:B,:A) => -1) #B => A
))
set_props!(net, 4, Dict(
    :initial_inventory => Dict(:A => 105), #initial inventory at storage
    :inventory_capacity => Dict(:A => Inf, :B => 0, :C => 0),
    :demand_distribution => Dict(:A => Normal(5,1)), #demand
    :demand_period => Dict(:A => 1) #order every other day on average
)) 
set_props!(net, 5, Dict(
    :initial_inventory => Dict(:A => 60), #initial inventory at retail
    :inventory_capacity => Dict(:A => Inf, :B => 0, :C => 0),
    :demand_distribution => Dict(:A => Normal(10,1)), #demand
    :demand_period => Dict(:A => 1) #order daily
))

##specify lead times
set_prop!(net, 1, 2, :lead_time, Dict(:C => Poisson(8)))
set_prop!(net, 2, 3, :lead_time, Dict(:B => Poisson(6)))
set_prop!(net, 3, 4, :lead_time, Dict(:A => Poisson(4)))
set_prop!(net, 4, 5, :lead_time, Dict(:A => Poisson(2)))

##define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
policy_variable = :echelon_stock #variable tracked by policy
s = Dict((2,:C) => 215, (3,:B) => 190, (4,:A) => 165, (5,:A) => 60) #lower bound on inventory
S = Dict((2,:C) => 215, (3,:B) => 190, (4,:A) => 165, (5,:A) => 60) #base stock level

##create environment and run simulation with reorder policy
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_type, review_period, policy_variable)

##make plots
using DataFrames, StatsPlots
#inventory level
df = filter(i -> i.node > 1, env.echelon_stock)
fig1 = @df df plot(:period, :level, group={Node = :node, Mat = :material}, linetype=:steppost, legend = :bottomright,
                    xlabel="period", ylabel="echelon stock", ylim=(min(0, minimum(df.level)),maximum(df.level)))