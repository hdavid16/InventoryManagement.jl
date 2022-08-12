#2-echelon system with variable demand
using Distributions
using InventoryManagement
using StatsBase

#define network topology and materials
# Node 1 = distribution center with unlimited inventory of A
# Node 2 = retail with market demand for A
net = MetaDiGraph(path_digraph(2))
set_prop!(net, :materials, [:A])

#parameters
D = Normal(5,1) #market demand distribution; mean = 5, std = 1
LT = Normal(5,1) #lead time distribution; mean = 5, std = 1
function base_stock(D,LT,service_target)
    z = quantile(Normal(),service_target) #z-score to reach target
    LTD = mean(D)*mean(LT+1) #expected demand over lead time (+1 is for delay in demand based on sequence of events)
    SS = z*sqrt(mean(D)^2*std(LT)^2 + mean(LT+1)*std(D)^2) #safety stock required (+1 is for delay in demand based on sequence of events)
    return LTD, SS #base stock
end
service = 0.9 #service level target
LTD, SS = base_stock(D,LT,service)
BS = LTD + SS

#specify lead time
set_prop!(net, 1, 2, :lead_time, Dict(:A => LT)) 

#specify node parameters
set_props!(net, 1, Dict(:initial_inventory => Dict(:A => Inf)))
set_props!(net, 2, Dict(
    :initial_inventory => Dict(:A => BS),
    :demand_distribution => Dict(:A => D), 
    :demand_frequency => Dict(:A => 1), #daily demand
    :market_partial_fulfillment => Dict(:A => false)
))

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods, backlog = true);

#define reorder policy parameters
policy_type = :sS #(s, S) policy
review_period = 1 #review every day periods
s = Dict((2,:A) => BS) #lower bound on inventory
S = Dict((2,:A) => BS) #base stock level
metrics=(:service_levels,) #calculate service level
window=(15,100) #exclude first 15 days from metrics evaluation

#run simulation with reorder policy
simulate_policy!(env, s, S; policy_type, review_period, window, metrics)
# service level on arc (2,:market) ~ 90%

##make plots
using DataFramesMeta, StatsPlots
backlog = @rsubset(env.inventory, :location == 2, :type == Symbol("unfulfilled"))
position = @rsubset(env.inventory, :location == 2, :type == Symbol("position"))
level = @rsubset(env.inventory, :location == 2, :type == Symbol("level"))
# on_hand = @rsubset(env.inventory, :location == 2, :type == Symbol("on_hand"))
fig1 = plot(
    [0,100],[BS, BS], 
    linestyle = :dash, linewidth = 2,
    label="base stock",
    xlabel="period", ylabel="level",
    title="Node 2", xlim=(0,100),
    palette=:seaborn_bright6,
    legend=:outertopright
)
plot!(fig1, 
    [0,100], [SS, SS], 
    linestyle = :dash, linewidth = 2,
    label="safety stock"
)
@df backlog plot!(fig1,
    :period, :amount, linewidth = 2,
    linetype=:steppost, label="backlog", 
)
@df position plot!(fig1,
    :period, :amount, linewidth = 2,
    linetype=:steppost, label="position"
)
@df level plot!(fig1,
    :period, :amount, linewidth = 2,
    linetype=:steppost, label="level"
)
# @df on_hand plot!(fig1,
#     :period, :amount,
#     linetype=:steppost, label="on-hand"
# )