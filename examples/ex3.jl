#3-echelon supply chain with market demand at levels 1 and 2 of the chain
using Distributions
using InventoryManagement

##define network topology and materials involved:
# Node 1 = Raw material supplier
# Nodes 2 = Plant
# Node 3 = Retail (direct demand of A occurs here)
net = MetaDiGraph(3)
connect_nodes!(net,
    1 => 2, 
    2 => 2, #production self-loop
    2 => 3,
)
set_prop!(net, :materials, [:A, :B, :C])

##specify node parameters
set_props!(net, 1, Dict(
    :initial_inventory => Dict(:C => Inf),
    :inventory_capacity => Dict(:A => 0, :B => 0, :C => Inf)
))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:A => 105, :B => 190, :C => 215), #initial inventory at plant
    :inventory_capacity => Dict(:A => Inf, :B => Inf, :C => Inf),
    :bill_of_materials => Dict((:C,:B) => -1, (:B,:A) => -1), #C => B; B => A
    :demand_distribution => Dict(:A => Normal(5,1)), #demand
    :demand_frequency => Dict(:A => 1), #order every other day on average,
    :supplier_priority => Dict(:A => 2, :B => 2, :C => 1)
))
set_props!(net, 3, Dict(
    :initial_inventory => Dict(:A => 60), #initial inventory at retail
    :inventory_capacity => Dict(:A => Inf, :B => 0, :C => 0),
    :demand_distribution => Dict(:A => Normal(10,1)), #demand
    :demand_frequency => Dict(:A => 1) #order daily
))

##specify lead times
set_prop!(net, 1, 2, :lead_time, Dict(:C => Poisson(8)))
set_prop!(net, 2, 2, :lead_time, Dict(:B => Poisson(6), :A => Poisson(4))) #production times
set_prop!(net, 2, 3, :lead_time, Dict(:A => Poisson(2)))

##define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #continuous review
policy_variable = :echelon_stock #variable tracked by policy
s = Dict((2,:C) => 215, (2,:B) => 190, (2,:A) => 165, (3,:A) => 60) #lower bound on inventory
S = Dict((2,:C) => 215, (2,:B) => 190, (2,:A) => 165, (3,:A) => 60) #base stock level

##create environment and run simulation with reorder policy
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true, reallocate = true)
simulate_policy!(env, s, S; policy_type, review_period, policy_variable)

##make plots
using DataFramesMeta, StatsPlots
#inventory level
fig1 = @chain env.inventory begin
    @rsubset(:type == Symbol("echelon"))
    @rsubset(:location > 1)
    plot(
        _.period, _.amount, group=(Node = _.location, Mat = _.material), 
        linetype=:steppost, legend = :bottomright,
        xlabel="period", ylabel="echelon stock", 
        ylim=(min(0, minimum(_.amount)), maximum(_.amount))
    )
end