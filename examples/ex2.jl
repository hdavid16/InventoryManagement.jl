#2-echelon system with variable demand

using Distributions
using InventoryManagement

#define network topology and materials
# Node 1 = distribution center with unlimited inventory of A
# Node 2 = retail with market demand for A
net = MetaDiGraph(path_digraph(2))
set_prop!(net, :materials, [:A])

#specify node parameters
set_props!(net, 1, Dict(:initial_inventory => Dict(:A => Inf)))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:A => 100),
    :demand_distribution => Dict(:A => Normal(5,0.5)),
    :demand_period => Dict(:A => 2))
)
#specify lead time
set_prop!(net, 1, 2, :lead_time, Dict(:A => Poisson(5)))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true)

#define reorder policy parameters
policy_type = :rQ #(r, Q) policy
review_period = 25 #review every 25 periods
r = Dict((2,:A) => 20) #lower bound on inventory
Q = Dict((2,:A) => 80) #base stock level

#run simulation with reorder policy
simulate_policy!(env, r, Q; policy_type, review_period)

#make plots
using StatsPlots
#unfulfilled market demand
demand = filter(i -> i.arc[1] == 2, env.demand)
fig1 = plot(demand.period, demand.unfulfilled, linetype=:steppost, lab="backlog",
                    xlabel="period", ylabel="level", title="Node 2, Material A")
#add inventory position
inventory_on_hand = filter(i -> i.level < Inf, env.inventory_on_hand)
plot!(fig1, inventory_on_hand.period, inventory_on_hand.level, linetype=:steppost, lab = "on-hand inventory")
