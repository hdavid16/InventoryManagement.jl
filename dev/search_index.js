var documenterSearchIndex = {"docs":
[{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/#Environment","page":"API","title":"Environment","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.SupplyChainEnv\nInventoryManagement.connect_nodes!\nInventoryManagement.reset!\nInventoryManagement.is_terminated","category":"page"},{"location":"api/#InventoryManagement.SupplyChainEnv","page":"API","title":"InventoryManagement.SupplyChainEnv","text":"Supply Chain Simulation Environment Constructor\n\nFields:\n\nnetwork::MetaDiGraph: Supply chain network directed graph.\nmarkets::Vector: Vector of market nodes where demand occurs.\nproducers::Vector: Vector of producer nodes where material transformation occurs.\nechelons::Dict: Dictionary with Vector of nodes downstream of each node in the network (including that node).\nmaterials::Vector: Vector with the names of all materials in the system.\ninventory_on_hand::DataFrame: Timeseries with on hand inventories @ each node.\ninventory_level::DataFrame: Timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)\ninventory_pipeline::DataFrame: Timeseries with pipeline inventories on each arc.\ninventory_position::DataFrame: Timeseries with inventory positions @ each node (inventory level + placed replenishments).\nechelon_stock::DataFrame: Timeseries with echelon inventory position @ each node.\ndemand::DataFrame: Timeseries with replenishment orders and market demand placed on each arc.\norders::DataFrame: History of all orders received.\nopen_orders::DataFrame: Temporary table with outstanding orders.\nshipments::DataFrame: Temporary table with active shipments and time to arrival on each arc.\nprofit::DataFrame: Timeseries with profit @ each node.\nmetrics::DataFrame: Service metrics (service level and fill rate) for each supplier and material.\nreward::Float64: Final reward in the system (used for RL)\nperiod::Int: Current period in the simulation.\nnum_periods::Int: Number of periods in the simulation.\nnum_orders::Int: Number of orders in the system.\ndiscount::Float64: Time discount factor (i.e. interest rate).\noptions::Dict: Simulation options\nbacklog::Bool: Indicator if backlogging is allowed.\nreallocate::Bool: Indicator if unfulfilled requests should be reallocated to alternate suppliers.\nevaluate_profit::Bool: Indicator if the profit should be evaluated at each node.\ncapacitated_inventory::Bool: Indicator if inventory limits should be enforced.\nguaranteed_service::Bool: Indicator if simulation should force lost sales after service time expires.\nseed::Int: Random seed.\n\n\n\n\n\n","category":"type"},{"location":"api/#InventoryManagement.connect_nodes!","page":"API","title":"InventoryManagement.connect_nodes!","text":"connect_nodes!(net::MetaDiGraph, arcs...)\n\nConnect nodes identified by arc pairs (e.g., 1 => 2).\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.reset!","page":"API","title":"InventoryManagement.reset!","text":"reset!(env::SupplyChainEnv)\n\nReset a SupplyChainEnv (empty all logging dataframes and set simulation time to 0).\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.is_terminated","page":"API","title":"InventoryManagement.is_terminated","text":"is_terminated(env::SupplyChainEnv)\n\nCheck if a simulation has terminated (i.e., has reached the maximum number of periods).\n\n\n\n\n\n","category":"function"},{"location":"api/#Inventory-Policy","page":"API","title":"Inventory Policy","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.reorder_policy\nInventoryManagement.simulate_policy!","category":"page"},{"location":"api/#InventoryManagement.reorder_policy","page":"API","title":"InventoryManagement.reorder_policy","text":"reorder_policy(env::SupplyChainEnv, reorder_point::Dict, policy_param::Dict,\n                    policy_type::Union{Dict, Symbol} = :rQ, \n                    review_period::Union{Int, AbstractRange, Vector, Dict} = 1,\n                    min_order_qty::Union{Real, Dict} = 0,\n                    adjust_expected_consumption::Bool = true)\n\nApply an inventory policy to specify the replinishment orders for each material     throughout the SupplyChainEnv.\n\nArguments\n\nenv: inventory management environment\nreorder_point: the s or r parameter in each node for each material in the system. The keys are of the form (node, material).\npolicy_param: the S or Q parameter in each node for each material in the system. The keys are of the form (node, material).\npolicy_type: :rQ for an (r,Q) policy, or :sS for an (s,S) policy. If passing a Dict, the policy type should be specified for each node (keys).\nreview_period: number of periods between each inventory review (Default = 1 for continuous review.). If a AbstractRange or Vector is used, the review_period indicates which periods the review is performed on. If a Dict is used, the review period should be specified for each (node, material) Tuple (keys). The values of this Dict can be either Int, AbstractRange, or Vector. Any missing (node, material) key will be assigned a default value of 1.\nmin_order_qty: minimum order quantity (MOQ) at each supply node. If a Dict is passed, the MOQ should be specified for each (node, material) Tuple (keys). The values should be Real. Any missing key will be assigned a default value of 0.\nadjust_expected_consumption: should the system be assumed to be centralized? If true then the upstream nodes know how much each downstream node is going to request and adjust the stock state to account for this.\n\n\n\n\n\n","category":"function"},{"location":"api/#InventoryManagement.simulate_policy!","page":"API","title":"InventoryManagement.simulate_policy!","text":"simulate_policy!(env::SupplyChainEnv, args...)\n\nStep through a simulation using a specified reorder policy. args are the arguments that are passed to the reorder_policy function.\n\n\n\n\n\n","category":"function"},{"location":"api/#Reorder-Actions","page":"API","title":"Reorder Actions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.show_action","category":"page"},{"location":"api/#InventoryManagement.show_action","page":"API","title":"InventoryManagement.show_action","text":"(x::SupplyChainEnv, action::Vector{T} where T <: Real)\n\nConvert a replenishment order vector into a NamedArray indicating how much of each material (rows) is being requested on each arc (columns).\n\n\n\n\n\n","category":"function"},{"location":"api/#Spaces","page":"API","title":"Spaces","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"InventoryManagement.action_space\nInventoryManagement.state\nInventoryManagement.state_space","category":"page"},{"location":"events/#Sequence-of-Events","page":"Sequence of Events","title":"Sequence of Events","text":"","category":"section"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"The sequence of events in each period of the simulation is patterned after that of the News Vendor Problem:","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"Start period.\nPlace inventory replenishment orders at each node by traversing the supply network downstream (using topological sorting).","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"If backlog = true, the previous period's backlog is added to the replenishment order. \nThe supplier to the node placing the order will attempt to fill the replenishment order via its on-hand inventory if possible. If the supplier is a producer node and its on-hand inventory is insufficient, the supplier will then attempt to fulfill the order via material production (if there is sufficient production capacity and raw material inventory). \nIf reallocate = true, then any amount that cannot be satisfied is reallocated to the next supplier in the supplier priority list (the lowest priority supplier will reallocate back to the highest priority supplier). \nAccepted replenishment orders are immediately shipped with a lead time sampled from the specified distribution. For distributor nodes, the lead time is the in-transit (transportation) time between distributor nodes. For producer nodes, the lead time is the plant production time. \nIf the lead time is 0, the stock will be immediately available to the requesting node so that it can be used to fulfill downstream orders as they arrive (possible due to the topological sorting).","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"Receive inventory that has arrived at each node (after the lead time has transpired).\nMarket demand for each material occurs after tossing a weighted coin with the probability of demand occurring defined by the inverse of the demand_period (average number of periods between positive external demands).","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"For example, if demand_period = 2, there is a 1/2 = 50% chance of having positive demand, or once every 2 days on average.","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"Demand (including any backlog if backlog = true) is fulfilled up to available inventory at the market nodes.\nUnfulfilled demand is backlogged (if backlog = true).\nAccounts for each node are generated:","category":"page"},{"location":"events/","page":"Sequence of Events","title":"Sequence of Events","text":"Accounts payable: invoice for fulfilled replenishment orders (payable to suppliers), invoice for delivered replenishment orders (payable to third-party shipper), pipeline inventory holding cost for in-transit inventory (cost to requestor), on-hand inventory holding cost, penalties for unfulfilled demand (cost to supplier).\nAccounts receivable: sales for internal and external demand.","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Example-1","page":"Examples","title":"Example 1","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example has plant with unlimited raw material supply that converts :B to :A with a 1:1 stoichiometry. The plant sells both materials to a downstream retailer that has market demand for both materials. This system is modeled using 3 nodes:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Plant: Node 1 (stores :B) => Node 2 (stores :A)\nRetailer: Node 3 buys :B from Node 1 and :A from Node 2","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Demand and lead times are deterministic. A continuous review (s,S) policy is used. 100 periods are simulated.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"See code with system and policy parameters here.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: ) (Image: )","category":"page"},{"location":"examples/#Example-2","page":"Examples","title":"Example 2","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example has a distributor with unlimited inventory (Node 1) that sells :A to a retailer with market demand (Node 2).","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Demand and lead time is stochastic. A periodic review (r,Q) policy is used. 100 periods are simulated.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"See code with system and policy parameters here.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"examples/#Example-3","page":"Examples","title":"Example 3","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example has plant that converts :C to :B to :A with a 1:1 stoichiometry for each reaction. The plant acquires raw materials from a supplier upstream with unlimited supply of :C and sells :A to a retailer downstream. There is direct market demand of :A at both the retailer and the plant. Thus, the plant has both internal and external demand. This system is modeled using 5 nodes:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Supplier: Node 1 (unlimited supply :C)\nPlant: Node 2 (stores raw :C) => Node 3 (stores intermediate :B) => Node 4 (stores product :A and sells it to both the retailer and the market)\nRetailer: Node 5 buys :A from Node 4 and sells :A to the market.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Demand and lead times are stochastic. A continuous review (s,S) policy is used. 100 periods are simulated.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"See code with system and policy parameters here.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"(Image: )","category":"page"},{"location":"model/#Model-Assumptions","page":"Supply Chain Model","title":"Model Assumptions","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The following assumptions hold in the current implementation, but can be modified in future releases.","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Production lead times are independent of the amount being produced.\nTransportation costs are paid to a third party (not a node in the network).\nReplenishment orders are placed in topological order. This means that upstream nodes place orders first. This allows the following scenario: if the lead time is 0, the ordered inventory will immediately be available so that the node can use it to fulfill downstream orders.","category":"page"},{"location":"model/#Model-Limitations","page":"Supply Chain Model","title":"Model Limitations","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The following features are not currently supported:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Alternate bills of materials (see Model Inputs) for the same material are not currently supported. This is particularly relevant for chemical systems (e.g., there are two reactions that can be used to make the same product).\nCapacity limitations on shared feedstock inventory among producer nodes (e.g., shared inventory tanks) are not enforced in a straightforward way since producer nodes have dedicated raw material storage. Shared feedstock inventory can be modeled by having an upstream storage node with zero lead time to each of the producer nodes. Each of the producer nodes should use an (s,S) replenishment policy with s = 0, S = 0. When a production order is of size x is going to be placed in a period, the policy will assume the feedstock position at the producer node is going to derop to -x and will order x to bring the position up to 0. Since the lead time is 0 and orders are placed with topological sorting (see Sequence of Events), the inventory will be immediately sent to the producer node and be available for when the production order comes in. However, if there is not enough production capacity to process x, the excess will be left in the dedicated storage for that producer node, making it possible to violate the shared inventory capacity constraint.\nIf a producer can produce more than 1 material, it is possible for it to produce all materials it is capable of producing simultaneously (if there are enough raw materials). This occurs because the model does not account for resource constraints (e.g., single reactor can only do reaction 1 or reaction 2, but not both simultaneously). However, these can be enforced manually with the reorder actions. Potential fixes (requires changing the source code):\nDrop inventory capacities to 0 when the production equipment is occupied. Requires modeling each production unit as its own node.\nDevelop a production model (perhaps based on the Resource-Task Network paradigm)","category":"page"},{"location":"model/#Creating-a-Supply-Chain-Network","page":"Supply Chain Model","title":"Creating a Supply Chain Network","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The supply network topology must be mapped on a network graph using Graphs.jl. The system parameters are stored in the network's metadata using MetaGraphs.jl. The network can be generated by using the MetaDiGraph function and passing one of the following:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Number of nodes in the network, which can the be connected via the connect_nodes! function:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"net = MetaDiGraph(3)\nconnect_nodes!(net,\n  1 => 2,\n  2 => 3\n)","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"An adjacency matrix for the network:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"net = MetaDiGraph(\n  [0 1 0;\n   0 0 1;\n   0 0 0]\n)","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"An existing DiGraph, such as a serial directed graph (path_digraph):","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"net = MetaDiGraph(path_digraph(3))","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"General network, node-specific, and arc-specific metadata can be added using the set_prop! and set_props! functions. The following subsections describe the property keys accepted for the supply chain network metadata.","category":"page"},{"location":"model/#General-Network-Parameters","page":"Supply Chain Model","title":"General Network Parameters","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The graph metadata should have the following fields in its metadata:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":":materials::Vector with a list of all materials in the system.","category":"page"},{"location":"model/#Node-specific-Parameters","page":"Supply Chain Model","title":"Node-specific Parameters","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Producers will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":":initial_inventory::Dict: initial inventory for each material (keys). Default = 0.\n:inventory_capacity::Dict: maximum inventory for each material (keys). Default = Inf.\n:holding_cost::Dict: unit holding cost for each material (keys). Default = 0.\n:early_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled before their due date for each material (keys). Default = true.\n:partial_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled partially for each material (keys). Default = true.\n:supplier_priority::Dict: (only when the node has at least 1 supplier) Vector of suppliers (from high to low priority) for each material (keys). When a request cannot be fulfilled due to insufficient production capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if reallocate = true). Default = inneighbors(SupplyChainEnv.network, node).\n:production_capacity::Dict: maximum production capacity for each material (keys). Default = Inf.\n:bill_of_materials::Union{Dict,NamedArray}: keys are material Tuples, where the first element is the input material and the second element is the product/output material; the values indicate the amount of input material consumed to produce 1 unit of output material. Alternatively, a NamedArray can be passed where the input materials are the rows and the output materials are the columns. The following convention is used for the bill of material (BOM) values:\nzero: input not involved in production of output.\nnegative number: input is consumed in the production of output.\npositive number: input is a co-product of the output.","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Distributors will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":":initial_inventory::Dict: initial inventory for each material (keys). Default = 0.\n:inventory_capacity::Dict: maximum inventory for each material (keys). Default = Inf.\n:holding_cost::Dict: unit holding cost for each material (keys). Default = 0.\n:early_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled before their due date for each material (keys). Default = true.\n:partial_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled partially for each material (keys). Default = true.\n:supplier_priority::Dict: (only when the node has at least 1 supplier) Vector of supplier priorities (from high to low) for each material (keys). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if reallocate = true). Default = inneighbors(SupplyChainEnv.network, node).","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Markets will have the following fields in their node metadata:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":":initial_inventory::Dict: initial inventory for each material (keys). Default = 0.\n:inventory_capacity::Dict: maximum inventory for each material (keys). Default = Inf.\n:holding_cost::Dict: unit holding cost for each material (keys). Default = 0.\n:early_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled before their due date for each material (keys). Default = true.\n:partial_fulfillment::Dict: (only when the node has at least 1 supplier) (true/false) on if the node accepts orders being fulfilled partially for each material (keys). Default = true.\n:supplier_priority::Dict: (only when the node has at least 1 supplier) Vector of supplier priorities (from high to low) for each material (keys). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if reallocate = true). Default = inneighbors(SupplyChainEnv.network, node).\n:demand_distribution::Dict: probability distributions from Distributions.jl for the market demands for each material (keys). For deterministic demand, instead of using a probability distribution, use D where D <: Number. Default = 0.\n:demand_period::Dict: mean number of periods between demand arrivals for each material (keys). Default = 1.\n:demand_sequence::Dict: a user specified Vector of market demand for each material (keys). When a nonzero Vector is provided, the demand_distribution and demand_period parameters are ignored. Default = zeros(SupplyChainEnv.num_periods).\n:sales_price::Dict: market sales price for each material (keys). Default = 0.\n:unfulfilled_penalty::Dict: unit penalty for unsatisfied market demand for each material (keys). Default = 0.\n:service_time::Dict: service time (probability distribution or deterministic value) allowed to fulfill market demand for each material (keys). Default = 0.","category":"page"},{"location":"model/#Arc-specific-Parameters","page":"Supply Chain Model","title":"Arc-specific Parameters","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"All arcs have the following fields in their metadata:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":":sales_price::Dict: unit sales price for inventory sent on that edge (from supplier to receiver) for each material (keys). Default = 0.\n:transportation_cost::Dict: unit transportation cost for shipped inventory for each material (keys). Default = 0.\n:pipeline_holding_cost::Dict: unit holding cost per period for inventory in-transit for each material (keys). Default = 0.\n:unfulfilled_penalty::Dict: unit penalty for unsatisfied internal demand for each material (keys). Default = 0.\n:lead_time::Distribution{Univariate, Discrete}: probability distributions from Distributions.jl for the lead times for each material (keys) on that edge. Lead times are transportation times when the edge has two distributor nodes and production times when the edge joins the producer and distributor nodes in a plant. For deterministic lead times, instead of using a probability distribution, use L where L  <: Number. Default = 0.\n:service_time::Dict: service time (probability distribution or deterministic value) allowed to fulfill internal demand for each material (keys). Default = 0.","category":"page"},{"location":"model/#Creating-a-Supply-Chain-Environment","page":"Supply Chain Model","title":"Creating a Supply Chain Environment","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The SupplyChainEnv function can be used to create a SupplyChainEnv Constructor. This function takes the following inputs:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"Positional Arguments:\nNetwork::MetaDiGraph: supply chain network with embedded metadata\nnum_periods::Int: number of periods to simulate\nKeyword Arguments (system options):\nbacklog::Bool = true: backlogging allowed if true; otherwise, orders that reach their due date and are not fulfilled become lost sales.\nreallocate::Bool = false: the system try to reallocate requests if they cannot be satisfied if true; otherwise, no reallocation is attempted.\nguaranteed_service::Bool = false: the simulation will operate under the assumptions in the Guaranteed Service Model (GSM). If true, backlog = true will be forced. Orders that are open and within the service time window will be backlogged. Once the service time expires, the orders become lost sales. In order to replicate the GSM assumption that extraordinary measures will be used to fulfill any expired orders, a dummy node with infinite supply can be attached to each node and set as the lowest priority supplier to that node.\ncapacitated_inventory::Bool = true: the simulation will enforce inventory capacity constraints by discarding excess inventory at the end of each period if true; otherwise, the system will allow the inventory to exceed the specified capacity.\nevaluate_profit::Bool = true: the simulation will calculate the proft at each node if true and save the results in SupplyChainEnv.profit.\nAditional Keyword Arguments:\ndiscount::Float64 = 0.: discount factor (i.e., interest rate) to account for the time-value of money.\nseed::Int = 0: random seed for simulation.","category":"page"},{"location":"model/#Simulation-Outputs","page":"Supply Chain Model","title":"Simulation Outputs","text":"","category":"section"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"The SupplyChainEnv Constructor has the following fields to store the simulation results in DataFrames:","category":"page"},{"location":"model/","page":"Supply Chain Model","title":"Supply Chain Model","text":"inventory_on_hand: on-hand inventory for each node, material, and period. Discarded inventory is marked when capacitated_inventory = true\ninventory_level: inventory level for each node, material, and period\ninventory_pipeline: in-transit inventory for each arc, material, and period\ninventory_position: inventory position for each node, material, and period\nechelon_stock: inventory position for each echelon, material, and period\ndemand: internal and external demands for each material on each arc (for internal demand) and each node (for external demand), at each period. The total demand quantities, fulfilled demand quantities, lead times and unfulfilled demand quantities are tabulated. If reallocate = true and the unfulfilled demand is reallocated, the arc that the demand is reallocated to is also indicated.\norders: internal and external orders for each material on each arc (for internal demand) and each node (for external demand). The ID, creation date, and quantity are stored for each order. The fulfilled column has a vector of Tuples that indicate the fulfillment time, supplier, and amount fulfilled. More than one Tuple will be shown if the order has been split.\nopen_orders: open (not yet fulfilled) internal and external orders for each material on each arc (for internal demand) and each node (for external demand). The ID, creation date, and quantity are stored for each open order. The due column indicates the time left until the order is due (as specified by the service_time).\nshipments: current in-transit inventory for each arc and material with remaining lead time.\nprofit: time-discounted profit for each node at each period.\nmetrics: service metrics (service level and fillrate) for each supplier and material.","category":"page"},{"location":"#InventoryManagement.jl:","page":"Home","title":"InventoryManagement.jl:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Discrete-time simulation environment for Inventory Management in Supply Chain Networks.","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InventoryManagement.jl allows modeling a multi-period multi-product supply chain network under stochastic stationary demand and stochastic lead times. A supply network can be constructed using the following types of nodes:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Producers: Nodes where inventory transformation takes place (e.g., raw materials are converted to intermediates or finished goods). Material transformation, including reactive systems with co-products, are modeled using Bills of Materials (see Model Inputs section).\nDistributors: Nodes where inventory is stored and distributed (e.g., distribution centers).\nMarkets: Nodes where end-customers place final product orders (e.g., retailer). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"These types of nodes can be used to model the following components of a supply network:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Manufacturing Plants: Plants are modeled using at least two nodes joined by a directed arc:\nRaw material storage node: Upstream producer node that stores raw materials and the plant's bill of materials.\nProduct storage node: Downstream distributor node that stores the materials produced at the plant.\nArc: the time elapsed between the consumption of raw materials and the production of goods (production time) is modeled with the arc lead time.\nDistribution Centers: DCs are modeled using distributor nodes.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Note: Any node can be marked as a market node to indicate that there is external demand for one or more materials stored at that node. This allows external demand at distribution centers or at manufacturing plants (either for raw materials or products).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The simplest network that can be modeled is one with a single retailer (distributor) with external demand (market) that is supplied by a warehouse (distributor). However, more complex systems can be modelled as well.","category":"page"},{"location":"","page":"Home","title":"Home","text":"When defining a supply network, a SupplyChainEnv object is created based on system inputs and network structure. This object can then be used to simulate the inventory dynamics under a stochastic environment and a specified inventory management policy. Stochasticity is modeled in the external demand quantities at the market nodes and in the lead times between connected nodes in the network. However, deterministic values can also be used for external demand or lead times if desired. In each period of the simulation, a decision-maker can specify inventory replenishnment orders throughout the network (refered to as actions), which consist of the inventory quantities requested by each node to each immediate supplier for each material in the system. If no action is taken during the simulation, the inventory levels will eventually be depleted by the external demand. Depending on the system configuration, unfulfilled external or internal demand can be either backlogged or considered a lost sale.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The SupplyChainEnv can also potentially be used in conjunction with ReinforcementLearning.jl to train a Reinforcement Learning agent that places replenishment orders (actions) throughout the network.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package generalizes and extends and the inventory management environment available in OR-Gym.","category":"page"},{"location":"#Dependencies","page":"Home","title":"Dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"InventoryManagement.jl relies primarily on the following packages:","category":"page"},{"location":"","page":"Home","title":"Home","text":"DataFrames.jl: Tabulate results.\nDistributions.jl: Define probability distributions for the lead times on the network arcs and the demand quantities at the market nodes.\nGraphs.jl: Define supply network topology\nMetaGraphs.jl: Specify system parameters for each node or arc, or for the system as a whole.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package can be installed with the Julia package manager. From the Julia REPL, type ] to enter the Pkg REPL mode and run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add InventoryManagement","category":"page"},{"location":"","page":"Home","title":"Home","text":"For the master branch, run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add https://github.com/hdavid16/InventoryManagement.jl","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Author: Hector D. Perez\nPosition: Ph. D. Candidate @ Carnegie Mellon University\nEmail: hdperez@cmu.edu\nYear: 2021","category":"page"},{"location":"policies/#Inventory-replenishment-policies","page":"Inventory Policies","title":"Inventory replenishment policies","text":"","category":"section"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"At each iteration in the simulation, an action can be provided to the system, which consists of the replenishment orders placed on every link in the supply network. This action must be of type Vector{Real} and must be nonnegative of the form: [Arc_1_Material_1, Arc_1_Material_2, ..., Arc_1Material_M, Arc_2_Material_1, ..., Arc_2_Material_M, ..., Arc_A_Material_1, ..., Arc_A_Material_M], where the ordering in the arcs is given by edges(SupplyChainEnv.network) and the ordering in the materials by SupplyChainEnv.materials.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"An action vector can be visualized as a NamedArray using show_action(SupplyChainEnv, action):","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"material ╲ arc │ :Arc_1  :Arc_2 ... :Arc_A\n───────────────┼──────────────────────────\n:Material_1    │  \n:Material_2    │  \n...            │\n:Material_M    │  ","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"The function reorder_policy can be used to implement an inventory reorder policy at each node based its inventory position or echelon stock. Reorder quantities are placed to the node's priority supplier. The reorder policy is applied for each material at each node in reverse topological order. This allows upstream nodes to determine their reorder quantities with information about the reorder quantities placed by their successors (relevant for producer nodes to ensure that raw material replenishments are synced with production orders). The two most common policies used in industry are the (s,S) and (r,Q) policies.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"The reorder_policy takes the following inputs and returns an action vector.","category":"page"},{"location":"policies/","page":"Inventory Policies","title":"Inventory Policies","text":"env::SupplyChainEnv: inventory management environment\nreorder_point::Dict: the s or r parameter in each node for each material in the system. The keys are of the form (node, material).\npolicy_param::Dict: the S or Q parameter in each node for each material in the system. The keys are of the form (node, material).\npolicy_type::Union{Symbol, Dict}: :rQ for an (r,Q) policy, or :sS for an (s,S) policy. If passing a Dict, the policy type should be specified for each node (keys).\nreview_period::Union{Int, AbstractRange, Vector, Dict}: number of periods between each inventory review (Default = 1 for continuous review.). If a AbstractRange or Vector is used, the review_period indicates which periods the review is performed on. If a Dict is used, the review period should be specified for each (node, material) Tuple (keys). The values of this Dict can be either Int, AbstractRange, or Vector. Any missing (node, material) key will be assigned a default value of 1.\nmin_order_qty::Union{Real, Dict}: minimum order quantity (MOQ) at each supply node. If a Dict is passed, the MOQ should be specified for each (node, material) Tuple (keys). The values should be Real. Any missing key will be assigned a default value of 0. \nadjust_expected_consumption::Bool: indicator if the reorder point should be increased (temporarilly) at a producer node by the expected raw material consumption for an expected incoming production order.","category":"page"}]
}
