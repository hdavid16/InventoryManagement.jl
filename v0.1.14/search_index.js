var documenterSearchIndex = {"docs":
[{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/#Environment","page":"API","title":"Environment","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.SupplyChainEnv\nInventoryManagement.reset!\nInventoryManagement.is_terminated","category":"page"},{"location":"api/#InventoryManagement.SupplyChainEnv","page":"API","title":"InventoryManagement.SupplyChainEnv","text":"Supply Chain Simulation Environment Constructor\n\nFields:\n\nnetwork::MetaDiGraph: Supply chain network directed graph.\nmarkets::Vector: Vector of market nodes where demand occurs.\nproducers::Vector: Vector of producer nodes where material transformation occurs.\ndistributors::Vector: Vector of distribution centers (excludes market nodes).\nmaterials::Vector: Vector with the names of all materials in the system.\nbill_of_materials::Array: Square matrix with BOM (rows = inputs, cols = ouputs); indices follow materials list; row with positive value is a co-product, row with negative value is a input (reactant).\ninv_on_hand::DataFrame: Timeseries with on hand inventories @ each node.\ninv_level::DataFrame: Timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)\ninv_pipeline::DataFrame: Timeseries with pipeline inventories on each arc.\ninv_position::DataFrame: Timeseries with inventory positions @ each node (inventory level + placed replenishments).\nreplenishments::DataFrame: Timeseries with replenishment orders placed on each arc.\nshipments::DataFrame: Temp table with active shipments and time to arrival on each arc.\nproduction::DataFrame: Temp table with active material production commited to an arc and time to ship.\ndemand::DataFrame: Timeseries with realization of demand at each market, and amounts sold, unfulfilled demand, and backlog.\nprofit::DataFrame: Timeseries with profit @ each node.\nreward::Float64: Final reward in the system (used for RL)\nperiod::Int: Current period in the simulation.\nnum_periods::Int: Number of periods in the simulation.\ndiscount::Float64: Time discount factor (i.e. interest rate).\nbacklog::Bool: Indicator if backlogging is allowed.\nreallocate::Bool: Indicator if unfulfilled requests should be reallocated to alternate suppliers.\nseed::Int: Random seed.\n\n\n\n\n\n","category":"type"},{"location":"api/#InventoryManagement.reset!","page":"API","title":"InventoryManagement.reset!","text":"reset!(env::SupplyChainEnv)\n\nReset a SupplyChainEnv (empty all logging dataframes and set simulation time to 0).\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.is_terminated","page":"API","title":"InventoryManagement.is_terminated","text":"is_terminated(env::SupplyChainEnv)\n\nCheck if a simulation has terminated (i.e., has reached the maximum number of periods).\n\n\n\n\n\n","category":"function"},{"location":"api/#Inventory-Policy","page":"API","title":"Inventory Policy","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.reorder_policy\nInventoryManagement.simulate_policy!","category":"page"},{"location":"api/#InventoryManagement.reorder_policy","page":"API","title":"InventoryManagement.reorder_policy","text":"reorder_policy(env::SupplyChainEnv, param1::Dict, param2::Dict,\n                    kind::Symbol = :rQ, review_period::Int = 1)\n\nApply an inventory policy to specify the replinishment orders for each material     throughout the SupplyChainEnv.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.simulate_policy!","page":"API","title":"InventoryManagement.simulate_policy!","text":"simulate_policy!(env::SupplyChainEnv, args...)\n\nStep through a simulation using a specified reorder policy. args are the arguments that are passed to the reorder_policy function.\n\n\n\n\n\n","category":"function"},{"location":"api/#Reorder-Actions","page":"API","title":"Reorder Actions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.place_requests!\nInventoryManagement.update_production!\nInventoryManagement.update_shipments!\nInventoryManagement.enforce_inventory_limits!\nInventoryManagement.update_position!\nInventoryManagement.simulate_markets!\nInventoryManagement.calculate_profit!","category":"page"},{"location":"api/#InventoryManagement.place_requests!","page":"API","title":"InventoryManagement.place_requests!","text":"place_requests!(x::SupplyChainEnv, act::Array, arcs::Vector)\n\nPlace inventory replenishment requests throughout the network.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.update_production!","page":"API","title":"InventoryManagement.update_production!","text":"update_production!(x::SupplyChainEnv)\n\nUpdate completed production at each producer node and send to the arc that it     was commited to.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.update_shipments!","page":"API","title":"InventoryManagement.update_shipments!","text":"update_shipments!(x::SupplyChainEnv)\n\nUpdate inventories throughout the network for arrived shipments.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.enforce_inventory_limits!","page":"API","title":"InventoryManagement.enforce_inventory_limits!","text":"enforce_inventory_limits!(x::SupplyChainEnv)\n\nDiscard any excess inventory (exceeding the inventory capacity at each node).\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.update_position!","page":"API","title":"InventoryManagement.update_position!","text":"update_positions!(x::SupplyChainEnv, n::Int, p::Any)\n\nUpdate inventory position and inventory level for material p at node n.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.simulate_markets!","page":"API","title":"InventoryManagement.simulate_markets!","text":"simulate_markets!(x::SupplyChainEnv)\n\nOpen markets, apply material demands, and update inventory positions.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.calculate_profit!","page":"API","title":"InventoryManagement.calculate_profit!","text":"calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)\n\nCalculate profit at each node in the network     (sales - production cost - purchase costs - transportation costs - holding costs).\n\n\n\n\n\n","category":"function"},{"location":"api/#Utilities","page":"API","title":"Utilities","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.show_action\nInventoryManagement.check_inputs","category":"page"},{"location":"api/#InventoryManagement.show_action","page":"API","title":"InventoryManagement.show_action","text":"show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)\n\nVisualize a replenishment order vector as a DataFrame indicating how much of which material is being requested on each arc.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.check_inputs","page":"API","title":"InventoryManagement.check_inputs","text":"check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,\n                mrkts::Vector, plants::Vector, mats::Vector, bom::Array,\n                num_periods::Int)\n\nCheck inputs when creating a SupplyChainEnv.\n\n\n\n\n\n","category":"function"},{"location":"events/#Sequence-of-Events","page":"Sequence of Events","title":"Sequence of Events","text":"","category":"section"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"The following sequence of events occurs in each period of the simulation:","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"Start period.\nPlace inventory replenishment orders at each node. These are limited to the available production capacity (supplier is a producer) or available inventory. If reallocate = true, then any amount that cannot be satisfied is reallocated to the next priority supplier.\nDistributors ship inventory.\nProducers attempt to satisfy the requested amount in the following order: 1) using any on-hand inventory, 2) using any uncommitted scheduled production (only relevant when there is co-production), and 3) manufacturing material. Materials scheduled for production have a production lead time. Production costs are incurred at the start of production.\nProducers send orders that have completed (after the production lead time).\nReceive inventory that has arrived at each node (after the lead time has transpired).\nPay suppliers for inventory received.\nPay shipper for inventory shipped.\nMarket demand occurs after tossing a weighted coin with the probability of demand occurring defined by the demand_frequency.\nDemand (including any backlog if backlog = true) is fulfilled up to available inventory at the markets.\nUnfulfilled demand is penalized and backlogged (if backlog = true).\nEach node pays a holding cost and a transportation cost for on-hand inventory and in-transit inventory at each period.","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Example-1","page":"Examples","title":"Example 1","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example is for a 100 period simulation of a supply network with one plant (node 1) that supplies a retailer (node 3), with stochastic demand for product :A. Node 3, has an alternate supplier, which is a distribution center (node 2). Node 3 prefers replenishing from the plant, which has a lower lead time. A (s,S) reorder policy is used at the retailer. When the on-hand level for material :A is depleted at the plant, the plant begins transforming raw material :B into :A. There is limited raw material supply at the plant. When the raw material stocks-out, node 3 switches to node 2 for its supply.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"See code here.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"examples/#Example-2","page":"Examples","title":"Example 2","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example is for a 100 period simulation of a supply network with one warehouse (node 1) that supplies a retailer (node 2) with stochastic demand for product :A. A (r,Q) reorder policy is used at the retailer every 25 periods.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"See code here.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"model/#Model","page":"Simulation Model","title":"Model","text":"","category":"section"},{"location":"model/#Model-Assumptions","page":"Simulation Model","title":"Model Assumptions","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"The following assumptions hold in the current implementation, but can be modified in future releases.","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"Producers produce material on demand (make-to-order policy).\nProducers can hold inventory. Downstream replenishment orders are fulfilled first with any on-hand inventory, and then via production only after there is no on-hand inventory left.\nReplenishment orders can only be satisfied with current on-hand inventory or available production capacity.\nCommited production orders count towards the inventory position of the downstream node, even if they haven't yet shipped (due to production lead time).\nProduction lead times are fixed and independent of the amount being produced.\nBacklogging is allowed at all nodes. When reallocation = true, the unfulfilled quantity from the lowest priority supplier is added to the highest priority supplier request in the next period. Note: Unfulfilled requests are currently penalized at market nodes only. However, this can be changed in a new release by passing an unfulfilled penalty parameter in th metadata of each node.\nTransportation costs are paid to a third party (not a node in the network).","category":"page"},{"location":"model/#Model-Limitations","page":"Simulation Model","title":"Model Limitations","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"The following features are not currently supported:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"Producers do not operate under a make-to-stock policy since any material produced gets shipped downstream. However, this can be accomodated by adding a dumby node downstream of the producer that holds inventory produced by the producer with zero lead time in between the nodes. Thus, using a proper reorder policy, the producer can act as a make-to-stock system that pushes inventory to the inventory holding node.\nProducers do not have market demand. However, this can be modelled by adding a market node downstream of the producer with zero lead time in between the nodes.\nAlternate bills of materials (see Model Inputs) for the same material are not currently supported. This is particularly relevant for chemical systems. However, the following workarounds can be done:\nIf the alternate reaction pathway has a byproduct, then the main product can be included as a co-product in the bill of materials of the byproduct. For example: A system with 5 materials (:A - :E) can have two ways to produce :A, :B + :C -> :A and :D -> :A + :E. The column for material :A can have the bill of material: [0 -1 -1 0 0]. The column for material :E can have the bill of materials: [1 0 0 -1 0]. However, :A will only be produced by the second pathway if a request for :E is made.\nMake a copy of the material to specify an alternate pathway. This will require specifying parameters for the copied material throughout the network.\nCapacity limitations on shared feedstock inventory among producer nodes (e.g., shared inventory tanks) cannot be enforced directly. This is because the shared inventory is its own node and feeds the inventory holding area in the producer node. Thus the total inventory is the inventory at the inventory node plus the inventory positions at the producers. Capacity limitations must be enforced manually via the reorder actions. Potential fixes: (requires changing the package code)\nMake the inventory capacity dynamic (a function of the producer inventory holding sites).\nSplit production site from feedstock inventory into two nodes. Requires updating model logic and behavior.\nIf a producer can produce more than 1 material, it is possible to produce all materials it is capable of producing simultaneously. This does not account for resource constraints (e.g., single reactor can only do reaction 1 or reaction 2, but not both simultaneously). However, these can be enforced manually with the reorder actions. Potential fixes: (requires changing the package code)\nDrop inventory capacities to 0 when the production equipment is occupied. Requires modeling each production unit as its own node.\nDevelop a production model (perhaps based on the Resource-Task Network paradigm)","category":"page"},{"location":"model/#Model-Inputs","page":"Simulation Model","title":"Model Inputs","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"The supply network topology must be mapped on a network graph using MetaGraphs.jl. The system parameters are stored in the network's metadata.","category":"page"},{"location":"model/#Node-specific","page":"Simulation Model","title":"Node-specific","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"Producers will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":":initial_inventory::Dict: initial inventory for each material (keys)\n:inventory_capacity::Dict: maximum inventory for each material (keys)\n:holding_cost::Dict: unit holding cost for each material (keys)\n:supplier_priority::Dict: (only when the node has at least 1 supplier) Vector of supplier priorities (from high to low) for each material (keys). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if env.reallocate == true).\n:production_cost::Dict: unit production cost for each material (keys)\n:production_capacity::Dict: maximum production capacity for each material (keys).\n:production_time::Dict: production lead time for each material (keys).","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"Distributors will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":":initial_inventory::Dict: initial inventory for each material (keys)\n:inventory_capacity::Dict: maximum inventory for each material (keys)\n:holding_cost::Dict: unit holding cost for each material (keys)\n:supplier_priority::Dict: Vector of supplier priorities (from high to low) for each material (keys). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if env.reallocate == true).","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"Markets will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":":initial_inventory::Dict: initial inventory for each material (keys)\n:inventory_capacity::Dict: maximum inventory for each material (keys)\n:holding_cost::Dict: unit holding cost for each material (keys)\n:supplier_priority::Dict: Vector of supplier priorities (from high to low) for each material (keys). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if env.reallocate == true).\n:demand_distribution::Dict: probability distributions from Distributions.jl for the market demands for each material (keys). For deterministic demand, instead of using a probability distribution, use [D] where D is a Number.\n:demand_frequency::Dict: probability that demand will occur (value between 0.0 and 1.0) for each material (keys)\n:demand_sequence::Dict: a user specified Vector of market demand for each material (keys). When a nonzero Vector is provided, the demand_distribution and demand_frequency parameters are ignored.\n:sales_price::Dict: market sales price for each material (keys)\n:demand_penalty::Dict: unit penalty for unsatisfied market demand for each material (keys)","category":"page"},{"location":"model/#Edge-specific","page":"Simulation Model","title":"Edge-specific","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"All edges have the following fields in their metadata:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":":sales_price::Dict: unit sales price for inventory sent on that edge (from supplier to receiver) for each material (keys)\n:transportation_cost::Dict: unit transportation cost per period for inventory in-transit for each material (keys)\n:lead_time::Distribution{Univariate, Discrete}: probability distributions from Distributions.jl for the lead times for each material (keys) on that edge. For deterministic lead times, instead of using a probability distribution, use [L] where L is a Number.","category":"page"},{"location":"model/#General-Network","page":"Simulation Model","title":"General Network","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"The graph metadata should have the following fields in its metadata:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":":materials::Vector with a list of all materials in the system.\n:bill_of_materials::Matrix: bill of materials indicating the production recipies for the materials in the system. The row numbers correspond to the input materials and the column numbers to the output materials. The numbering matches that of the materials vector. The magnitude of each element is proportional to the production of one unit of output material. Each element can have one of three types of values:\nzero: input not involved in production of output.\nnegative number: input is consumed in the production of output.\npositive number: input is a co-product of the output.","category":"page"},{"location":"model/#Model-Output","page":"Simulation Model","title":"Model Output","text":"","category":"section"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"A SupplyChainEnv has the following fields:","category":"page"},{"location":"model/","page":"Simulation Model","title":"Simulation Model","text":"network::MetaDiGraph: Supply Chain Network (metagraph)\nmarkets::Array: list of market nodes\nproducers::Array: list of producer nodes\ndistributors::Array: list of distribution nodes (excludes end distributors where markets exist)\nmaterials::Array: list of all material (material) names (strings)\nbill_of_materials::Matrix square matrix with BOM (rows = input materials, cols = output materials; indices follow materials list; positive value is a co-product, negative is a feedstock)\ninv_on_hand::DataFrame: timeseries with on hand inventories @ each node.\ninv_level::DataFrame: timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)\ninv_pipeline::DataFrame: timeseries with pipeline inventories on each arc.\ninv_position::DataFrame: timeseries with inventory positions @ each node (inventory level + placed replenishments).\nreplenishments::DataFrame: timeseries Replenishment orders placed on each edge at the end of each period\nshipments::DataFrame: current shipments and time to arrival for each node\nproduction::DataFrame: current material production committed to an edge and lead time to ship. Note: byproducts are scheduled to go to the producing node n (edge (n,n)).\ndemand::DataFrame: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market\nprofit::DataFrame: timeseries with profit at each node\nreward::Float64: reward in the system (used for RL)\nperiod::Int: period in the simulation\nnum_periods::Int: number of periods in the simulation\ndiscount::Float64: time discount factor (interest rate)\nbacklog::Bool: backlogging allowed if true; otherwise, unfulfilled demand is lost sales\nreallocate::Bool: the system try to reallocate requests if they cannot be satisfied if true; otherwise, no reallocation is attempted.\nseed::Int: random seed","category":"page"},{"location":"#InventoryManagement.jl:","page":"Home","title":"InventoryManagement.jl:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Discrete-time simulation environment for Inventory Management in Supply Networks.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InventoryManagement.jl allows modeling a multi-period multi-product supply network. A supply network can be constructed using the following node types:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Producers: Nodes where inventory transformation takes place (e.g., intermediates or final materials are produced). Reactive systems, including those with co-products, can be modelled using Bills of Materials (see Model Inputs).\nDistributors: Intermediate nodes where inventory is stored and distributed (e.g., distribution centers).\nMarkets: Nodes where end-customers place final product orders (i.e., retailer). These are the last (sink) nodes in the network.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The simplest network that can be modeled is one with a single market with one producer or distributor. However, more complex systems can be modelled as well.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When defining a supply network, a SupplyChainEnv object is created based on system inputs and network structure. This object can then be used to execute a simulation of the inventory dynamics. During a simulation, stochastic demand at each of the markets can occur for each of the materials in each period. When product demand occurs at the market, sales are made based on available inventory. Any unfulfilled demand is either backlogged or considered a lost sale depending on the system definition. If no action is taken duirng the simulation, the inventory levels will eventually be depleted. To avoid this from happening, a decision-maker can interact with the system in each period by making inventory replenishment decisions (refered to as actions). Lead times for in-transit inventory as well as production lead times are accounted for in the simulation. Transportation lead times can be modelled stochastically to account for lead time uncertainty. From a service time perspective, demand at market nodes has zero service time, whereas non-market nodes have service time equal to the production lead time + transportation lead time.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The SupplyChainEnv can also potentially be used in conjunction with ReinforcementLearning.jl to train a Reinforcement Learning agent.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package generalizes and extends and the inventory management environment available in OR-Gym.","category":"page"},{"location":"#Dependencies","page":"Home","title":"Dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InventoryManagement.jl mainly relies on the following packages:","category":"page"},{"location":"","page":"Home","title":"Home","text":"MetaGraphs.jl: Define supply network structure and specify node- and edge-specific parameters.\nDataFrames.jl: Tabulate results.\nDistributions.jl: Define probability distributions for the lead times in between nodes and the demands at the market nodes.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add InventoryManagement","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Author: Hector D. Perez\nPosition: Ph. D. Candidate @ Carnegie Mellon University\nEmail: hdperez@cmu.edu\nYear: 2021","category":"page"},{"location":"policies/#Inventory-replenishment-policies","page":"Inventory Policies","title":"Inventory replenishment policies","text":"","category":"section"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"At each iteration in the simulation, an action can be provided to the system, which consists of the replenishment orders placed on every link in the supply network. This action must be of type Vector{Real} and must be nonnegative of the form: [Edge1_Material1, Edge1_Material2, ..., Edge1_MaterialM, Edge2_Material1, ..., Edge2_MaterialM, ..., EdgeE_Material1, ..., EdgeE_MaterialM], where the ordering in the edges is given by edges(env.network) and the ordering in the materials by env.materials.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"An action vector can be visualized as a DataFrame using show_action(action, env::SupplyChainEnv).","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"The function reorder_policy can be used to implement an inventory reorder policy at each node based its inventory position. Reorder quantities are placed to the node's priority supplier. The two most common policies used in industry are the (s,S) and (r,Q) policies.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"The reorder_policy takes the following inputs and returns an action vector.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"env::SupplyChainEnv: inventory management environment\nparam1::Dict: the s or r parameter in each node for each material in the system. The keys are of the form (node, material).\nparam2::Dict: the S or Q parameter in each node for each material in the system. The keys are of the form (node, material).\nkind::Symbol: :rQ for an (r,Q) policy, or :sS for an (s,S) policy\nreview_period::Int: number of periods between each inventory review (Default = 1 for continuous review.)","category":"page"}]
}
