using LightGraphs, MetaGraphs, Distributions
using InventoryManagement, StatsPlots, Random

#define network connectivity
net = MetaDiGraph(11)
# 1 = Heater (includes storage for A)
# 2 = Tank Hot A
# 3 = Tank B
# 4 = Tank C
# 5 = Reactor 1
# 6 = Market for Product 1
# 7 = Reactor 2
# 8 = Tank AB
# 9 = Tank BC
# 10 = Still (includes storage for E)
# 11 = Market for Product 2
add_edge!(net, 1, 2)
rxfeeds = [2, 3, 4, 8, 9] #reactor feeds (including recycle)
rxprods = [6, 8, 9, 10] #reactor products
for m in rxfeeds
    add_edge!(net, m, 5)
    add_edge!(net, m, 7)
end
for m in rxprods
    add_edge!(net, 5, m)
    add_edge!(net, 7, m)
end
add_edge!(net, 10, 11)
add_edge!(net, 10, 8)

#define materials
materials = [:A, :HotA, :B, :C, :AB, :BC, :E, :P1, :P2]
set_prop!(net, :materials, materials)
#              -1.0*A -> 1.0*HotA (column 2; heating)
#    -0.5*B + -0.5*C  -> 1.0*BC (column 6; reaction 1)
#    -0.2*C + -0.8*AB -> 1.0*E (column 7; reaction 3)
# -1.0*HotA + -1.5*BC -> 1.5*AB + 1.0*P1 (column 8; reaciton 2)
#          -1.0/0.9*E -> 0.1/0.9*AB + 1.0*P2 (column 9; separation)
bom = [0  -1  0  0  0     0     0     0         0;
       0   0  0  0  0     0     0  -1.0         0;
       0   0  0  0  0  -0.5     0     0         0;
       0   0  0  0  0  -0.5  -0.2     0         0;
       0   0  0  0  0     0  -0.8  +1.5  +0.1/0.9;
       0   0  0  0  0     0     0  -1.5         0;
       0   0  0  0  0     0     0     0  -1.0/0.9;
       0   0  0  0  0     0     0     0         0;
       0   0  0  0  0     0     0     0         0]
set_prop!(net, :bill_of_materials, bom)

#unlimited feeds (:A, :B, :C) and final product (:P1, :P2) storage

#specify producers
#initialize properties; ignore costs and prices
param = Dict(:initial_inventory => Dict(m => 0. for m in materials),
             :inventory_capacity => Dict(m => 0. for m in materials),
             :holding_cost => Dict(m => 0 for m in materials),
             :production_cost => Dict(m => 0 for m in materials),
             :production_time => Dict(m => 0 for m in materials),
             :production_capacity => Dict(m => 0. for m in materials),
             :supplier_priority => Dict(m => [] for m in materials))
#Heater (holds material :A)
heat = deepcopy(param)
heat[:initial_inventory][:A] = Inf
heat[:inventory_capacity][:A] = Inf
heat[:production_capacity][:HotA] = 100
heat[:production_time][:HotA] = 1
delete!(heat, :supplier_priority)
set_props!(net, 1, heat)
#Reactors (don't hold any inventory)
rx = deepcopy(param)
rx[:inventory_capacity][:HotA] = 80*0.4
rx[:inventory_capacity][:BC] = 80*0.6
rx[:inventory_capacity][:AB] = 80*0.8
rx[:inventory_capacity][:B] = 80*0.5
rx[:inventory_capacity][:C] = 80*0.5
rx[:production_capacity][:BC] = 80
rx[:production_time][:BC] = 2
rx[:production_capacity][:E] = 80
rx[:production_time][:E] = 1
rx[:production_capacity][:AB] = 80*0.4
rx[:production_capacity][:P1] = 80*0.6
rx[:production_time][:P1] = 2
suppliers_rx = inneighbors(net, 5)
rx[:supplier_priority][:HotA] = [2]
rx[:supplier_priority][:B] = [3]
rx[:supplier_priority][:C] = [4]
rx[:supplier_priority][:AB] = [8]
rx[:supplier_priority][:BC] = [9]
set_props!(net, 5, rx)
rx[:inventory_capacity][:HotA] = 50*0.4
rx[:inventory_capacity][:BC] = 50*0.6
rx[:inventory_capacity][:AB] = 50*0.8
rx[:inventory_capacity][:B] = 50*0.5
rx[:inventory_capacity][:C] = 50*0.5
rx[:production_capacity][:BC] = 50
rx[:production_capacity][:E] = 50
rx[:production_capacity][:AB] = 50*0.4
rx[:production_capacity][:P1] = 50*0.6
set_props!(net, 7, rx)
#Still (holds material :E)
still = deepcopy(param)
# still[:initial_inventory][:E] = 100
still[:inventory_capacity][:E] = 100
still[:production_time][:P2] = 2
still[:production_capacity][:P2] = 200*0.9
still[:production_capacity][:AB] = 200*0.1
suppliers_still = inneighbors(net, 10)
for p in materials #priority is irrelevant
    still[:supplier_priority][p] = suppliers_still
end
set_props!(net, 10, still)

#specify feed and intermediate tanks; ignore costs and prices
param = Dict(:initial_inventory => Dict(m => 0. for m in materials),
             :inventory_capacity => Dict(m => 0. for m in materials),
             :holding_cost => Dict(m => 0 for m in materials),
             :supplier_priority => Dict(m => [] for m in materials))
#HotA
hota = deepcopy(param)
# hota[:initial_inventory][:HotA] = 100
hota[:inventory_capacity][:HotA] = 100
suppliers_hota = inneighbors(net, 2)
for p in materials #priority is irrelevant
    hota[:supplier_priority][p] = suppliers_hota
end
set_props!(net, 2, hota)
#B
b = deepcopy(param)
b[:initial_inventory][:B] = Inf
b[:inventory_capacity][:B] = Inf
set_props!(net, 3, b)
#C
c = deepcopy(param)
c[:initial_inventory][:C] = Inf
c[:inventory_capacity][:C] = Inf
set_props!(net, 4, c)
#AB
ab = deepcopy(param)
# ab[:initial_inventory][:AB] = 200
ab[:inventory_capacity][:AB] = 200
suppliers_ab = inneighbors(net, 8)
for p in materials #priority is irrelevant
    ab[:supplier_priority][p] = suppliers_ab
end
set_props!(net, 8, ab)
#BC
bc = deepcopy(param)
# bc[:initial_inventory][:BC] = 150
bc[:inventory_capacity][:BC] = 150
suppliers_bc = inneighbors(net, 9)
for p in materials #priority is irrelevant
    bc[:supplier_priority][p] = suppliers_bc
end
set_props!(net, 9, bc)

#specify product tanks (markets); ignore costs and prices
param = Dict(:initial_inventory => Dict(m => 0. for m in materials),
             :inventory_capacity => Dict(m => 0. for m in materials),
             :holding_cost => Dict(m => 0 for m in materials),
             :demand_distribution => Dict{Symbol,Any}(m => [0] for m in materials),
             :demand_frequency => Dict(m => 1/7 for m in materials), #demand occurs once/week on average
             :sales_price => Dict(m => 0 for m in materials),
             :demand_penalty => Dict(m => 0 for m in materials),
             :supplier_priority => Dict(m => [] for m in materials))
#P1
p1 = deepcopy(param)
# p1[:initial_inventory][:P1] = 50
p1[:inventory_capacity][:P1] = Inf
p1[:demand_distribution][:P1] = Normal(10, 5)
suppliers_p1 = inneighbors(net, 6)
for p in materials #priority is irrelevant
    p1[:supplier_priority][p] = suppliers_p1
end
set_props!(net, 6, p1)
#P2
p2 = deepcopy(param)
# p2[:initial_inventory][:P2] = 50
p2[:inventory_capacity][:P2] = Inf
p2[:demand_distribution][:P2] = Normal(10, 5)
suppliers_p2 = inneighbors(net, 11)
for p in materials #priority is irrelevant
    p2[:supplier_priority][p] = suppliers_p2
end
set_props!(net, 11, p2)

#specify links
param = Dict(:sales_price => Dict(m => 0 for m in materials),
             :transportation_cost => Dict(m => 0 for m in materials),
             :lead_time => [0])
for e in edges(net)
    set_props!(net, e, param)
end
# #Heater -> HotA
# h_ha = deepcopy(param)
# h_ha[:lead_time] = [1]
# set_props!(net, 1, 2, h_ha)
# #Reactor (1&2) -> P1
# r_p1 = deepcopy(param)
# r_p1[:lead_time] = [2]
# set_props!(net, 5, 6, r_p1)
# set_props!(net, 7, 6, r_p1)
# #Reactor (1&2) -> AB
# r_ab = deepcopy(param)
# r_ab[:lead_time] = [2]
# set_props!(net, 5, 8, r_ab)
# set_props!(net, 7, 8, r_ab)
# #Reactor (1&2) -> BC
# r_bc = deepcopy(param)
# r_bc[:lead_time] = [2]
# set_props!(net, 5, 9, r_bc)
# set_props!(net, 7, 9, r_bc)
# #Reactor (1&2) -> Still (E)
# r_e = deepcopy(param)
# r_e[:lead_time] = [1]
# set_props!(net, 5, 10, r_e)
# set_props!(net, 7, 10, r_e)
# #Still -> AB
# s_ab = deepcopy(param)
# s_ab[:lead_time] = [1]
# set_props!(net, 10, 8, s_ab)
# #save param to all other links
# for e in edges(net)
#     if !(e in keys(net.eprops))
#         set_props!(net, e, param)
#     end
# end

#create environment
num_periods = 100
env = SupplyChainEnv(net, num_periods)
Random.seed!(env) #set random seed

#define reorder policy parameters
policy = :sS #(s, S) policy
on = :position #monitor inventory position
#initialize
s = Dict((n,m) => -1. for n in vertices(net), m in materials) #default is no trigger
S = Dict((n,m) => 0. for n in vertices(net), m in materials) #default is no reorder
#Reactors
for p in [:HotA,:B,:C,:AB,:BC]
    max_inv_rx1 = get_prop(net, 5, :inventory_capacity)[p]
    max_inv_rx2 = get_prop(net, 7, :inventory_capacity)[p]
    s[5,p], S[5,p] = max_inv_rx1, max_inv_rx1 #maintain full levels to have reactor ready for reaction when needed
    s[7,p], S[7,p] = max_inv_rx2, max_inv_rx2 #maintain full levels to have reactor ready for reaction when needed
end
#Still
max_E_tnk = get_prop(net, 10, :inventory_capacity)[:E]
s[10,:E], S[10,:E] = max_E_tnk, max_E_tnk #maintain full levels to have reactor ready for reaction when needed
#Markets
for (p, n) in zip([:P1,:P2],[6,11])
    s[(n,p)] = 60 #trigger at 5x mean demand
    S[(n,p)] = 60 #keep levels at 5x mean demand
end
#run simulation with reorder policy (set policy levels dynamically for the tanks)
for t in 1:env.num_periods
    #Tanks
    for (p, n) in zip([:HotA,:AB,:BC],[2,8,9])
        max_inv_tnk = get_prop(net, n, :inventory_capacity)[p]
        inv_rx1 = filter(i -> i.period == env.period && i.material == p && i.node == 5, env.inv_position).level[1]
        inv_rx2 = filter(i -> i.period == env.period && i.material == p && i.node == 7, env.inv_position).level[1]
        inv_sep = filter(i -> i.period == env.period && i.material == p && i.node == 10, env.inv_position).level[1]
        max_tnk = max_inv_tnk - inv_rx1 - inv_rx2 - inv_sep #level at tank cannot exceed this value
        if n == 2
            s[n,p], S[n,p] = max_tnk, max_tnk #keep level of HotA full
        else
            max_inv_rx1 = get_prop(net, 5, :inventory_capacity)[p]
            max_inv_rx2 = get_prop(net, 7, :inventory_capacity)[p]
            max_inv_sep = get_prop(net, 10, :inventory_capacity)[p]
            levl = min(max_tnk, max_inv_rx1 + max_inv_rx2 + max_inv_sep) #keep enough to feed all units (limit by max_tnk)
            s[n,p], S[n,p] = levl, levl #set levels
        end
    end
    #set action
    action = reorder_policy(env, s, S, on, policy, :priority)
    #step
    (env)(action)
end

#make plots
using DataFrames
#intermediate tanks
on_hand = filter(i -> i.level < Inf && i.material in [:AB, :BC, :E, :HotA], env.inv_on_hand)
on_hand_gb = groupby(on_hand, [:material, :period])
inv_on_hand = combine(on_hand_gb, :level => sum)
fig1 = @df inv_on_hand plot(:period, :level_sum, group=:material, linetype=:steppost, legend = :topleft,
                          xlabel="period", ylabel="tank level", ylims = [0,150])

#product tanks
on_hand = filter(i -> i.level < Inf && i.material in [:P1, :P2], env.inv_on_hand)
on_hand_gb = groupby(on_hand, [:material, :period])
inv_on_hand = combine(on_hand_gb, :level => sum)
fig2 = @df inv_on_hand plot(:period, :level_sum, group=:material, linetype=:steppost, legend = :topleft,
                          xlabel="period", ylabel="tank level")
#sales
sales = filter(i -> i.material == :P1 ? i.node == 6 :
                    i.material == :P2 ? i.node == 11 : false, env.demand)
sales_grouped = groupby(sales, :material)
sales = transform(sales_grouped, :sale => cumsum)

fig3 = @df sales plot(:period, :sale_cumsum, group={CumSales = :material}, linetype=:steppost, legend = :topleft,
                        xlabel="period", ylabel="amount")
@df sales plot!(:period, :backlog, group={Backlog = :material}, linetype=:steppost)
