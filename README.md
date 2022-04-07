# InventoryManagement.jl:

*Discrete-time simulation environment for Inventory Management in Supply Networks.*

[![DOI](https://zenodo.org/badge/357413023.svg)](https://zenodo.org/badge/latestdoi/357413023)

![](logo.png)

## Table of Contents

1. [Overview](#overview)
2. [Dependencies](#dependencies)
3. [Installation](#installation)
4. [Sequence of Events](#sequence-of-events)
5. [Model Assumptions](#model-assumptions)
6. [Model Limitations](#model-limitations)
7. [Inventory replenishment policies](#inventory-replenishment-policies)
8. [Model Inputs](#model-inputs)
    - [Node-specific](#node-specific)
    - [Edge-specific](#edge-specific)
    - [General (Network)](#general-network)
9. [Model Output: SupplyChainEnv](#model-output)
10. [Examples](#examples)
    - [Example #1: alternate suppliers and continuous review (s,S) policy](#example-1)
    - [Example #2: unlimited upstream supply and periodic (r,Q) policy](#example-2)
    - [Example #3: chemical production facility](#example-3)
11. [Contact](#contact)

## Overview

*InventoryManagement.jl* allows modeling a multi-period multi-product supply network under stochastic stationary demand and stochastic lead times. A supply network can be constructed using the following types of nodes:
- `Producers`: Nodes where inventory transformation takes place (e.g., raw materials are converted to intermediates or finished goods). Material transformation, including reactive systems with co-products, are modeled using [Bills of Materials](https://en.wikipedia.org/wiki/Bill_of_materials) (see [Model Inputs section](#node-specific)).
- `Distributors`: Nodes where inventory is stored and distributed (e.g., distribution centers).
- `Markets`: Nodes where end-customers place final product orders (e.g., retailer). 

These types of nodes can be used to model the following components of a supply network:
- `Manufacturing Plants`: Plants are modeled using two nodes joined by a directed arc:
  - Raw material storage node: Upstream `producer` node that stores raw materials and the plant's `bill of materials`.
  - Product storage node: Downstream `distributor` node that stores the materials produced at the plant.
  - Arc: the time elapsed between the consumption of raw materials and the production of goods (production time) is modeled with the arc lead time.
- `Distribution Centers`: DCs are modeled using `distributor` nodes.

Note: Any node can be marked as a `market` node to indicate that there is external demand for one or more materials stored at that node. This allows external demand at distribution centers or at manufacturing plants (either for raw materials or products).

The simplest network that can be modeled is one with a single retailer (`distributor`) with external demand (`market`) that is supplied by a warehouse (`distributor`). However, more complex systems can be modelled as well.

When defining a supply network, a `SupplyChainEnv` object is created based on system inputs and network structure. This object can then be used to simulate the inventory dynamics under a stochastic environment and a specified inventory management policy. Stochasticity is modeled in the external demand quantities at the `market` nodes and in the lead times between connected nodes in the network. However, deterministic values can also be used for external demand or lead times if desired. In each period of the simulation, a decision-maker can specify inventory replenishnment orders throughout the network (refered to as `actions`), which consist of the inventory quantities requested by each node to each immediate supplier for each material in the system. If no action is taken during the simulation, the inventory levels will eventually be depleted by the external demand. Depending on the system configuration, unfulfilled external or internal demand can be either backlogged or considered a lost sale.

The `SupplyChainEnv` can also potentially be used in conjunction with [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl) to train a Reinforcement Learning `agent` that places replenishment orders (`actions`) throughout the network.

This package generalizes and extends and the inventory management environment available in [OR-Gym](https://github.com/hubbs5/or-gym).

## Dependencies

*InventoryManagement.jl* relies primarily on the following packages:
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl): Tabulate results.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Define probability distributions for the lead times on the network arcs and the demand quantities at the market nodes.
- [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl): Define supply network topology
- [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl): Specify system parameters for each node or arc, or for the system as a whole.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```julia
pkg> add InventoryManagement
```

For the master branch, run:
```julia
pkg> add https://github.com/hdavid16/InventoryManagement.jl
```

## Sequence of Events

The sequence of events in each period of the simulation is patterned after that of the [News Vendor Problem](https://optimization.cbe.cornell.edu/index.php?title=Newsvendor_problem):
1. Start period.
2. Place inventory replenishment orders at each node by traversing the supply network upstream (using reverse [topological sorting](https://en.wikipedia.org/wiki/Topological_sorting)). If `SupplyChainEnv.options[:backlog] = true`, the previous period's backlog is added to the replenishment order. The supplier to the node placing the order will attempt to fill the replenishment order via its on-hand inventory if possible. If the supplier is a `producer` node and its on-hand inventory is insufficient, the supplier will then attempt to fulfill the order via material production (if there is sufficient `production capacity` and `raw material inventory`). If `SupplyChainEnv.options[:reallocate] = true`, then any amount that cannot be satisfied is reallocated to the next supplier in the `supplier priority` list. Accepted replenishment orders are immediately shipped with a lead time sampled from the specified distribution. For `distributor` nodes, the lead time is the in-transit (transportation) time between `distributor` nodes. For `producer` nodes, the lead time is the plant production time.
4. Receive inventory that has arrived at each node (after the lead time has transpired).
5. Market demand for each material occurs after tossing a weighted coin with the probability of demand occurring defined by the inverse of the `demand_period` (average number of periods between positive external demands, meaning that if `demand_period = 2`, there is a `1/2 = 50%` chance of having positive demand, or once every 2 days on average).
8. Demand (including any backlog if `SupplyChainEnv.options[:backlog] = true`) is fulfilled up to available inventory at the `market` nodes.
9. Unfulfilled demand is backlogged (if `SupplyChainEnv.options[:backlog] = true`) and penalized.
10. Accounts for each node are generated:
  - Accounts payable: invoice for fulfilled replenishment orders (payable to suppliers), invoice for delivered replenishment orders (payable to third-party shipper), pipeline inventory holding cost for in-transit inventory (cost to requestor), on-hand inventory holding cost, penalties for unfulfilled demand (cost to supplier)
  - Accounts receivable: sales for internal and external demand

## Model Assumptions

The following assumptions hold in the current implementation, but can be modified in future releases.

- `Producers` produce material on demand ([make-to-order](https://en.wikipedia.org/wiki/Build_to_order) policy).
- `Producers` can hold inventory. Downstream replenishment orders are fulfilled first with any on-hand inventory, and then via production only after there is no on-hand inventory left.
- Replenishment orders can only be satisfied with current on-hand inventory or available production capacity.
- Commited production orders count towards the inventory position of the downstream node, even if they haven't yet shipped (due to production lead time).
- Production lead times are fixed and independent of the amount being produced.
- Backlogging is allowed at all nodes. When `reallocation = true`, the unfulfilled quantity from the lowest priority supplier is added to the highest priority supplier request in the next period. *Note:* Unfulfilled requests are currently penalized at `market` nodes only. However, this can be changed in a new release by passing an unfulfilled penalty parameter in th metadata of each node.
- Transportation costs are paid to a third party (not a node in the network).

## Model Limitations

The following features are not currently supported:

- `Producers` do not operate under a [make-to-stock](https://en.wikipedia.org/wiki/Build_to_stock) policy since any material produced gets shipped downstream. However, this can be accomodated by adding a dumby node downstream of the `producer` that holds inventory produced by the `producer` with zero lead time in between the nodes. Thus, using a proper reorder policy, the `producer` can act as a `make-to-stock` system that pushes inventory to the inventory holding node.
- `Producers` do not have market demand. However, this can be modelled by adding a `market` node downstream of the `producer` with zero lead time in between the nodes.
- Alternate bills of materials (see [Model Inputs](#graph-specific)) for the same material are not currently supported. This is particularly relevant for chemical systems. However, the following workarounds can be done:
  - If the alternate reaction pathway has a byproduct, then the main product can be included as a co-product in the bill of materials of the byproduct. For example: A system with 5 materials (`:A - :E`) can have two ways to produce `:A`, `:B + :C -> :A` and `:D -> :A + :E`. The column for material `:A` can have the bill of material: `[0 -1 -1 0 0]`. The column for material `:E` can have the bill of materials: `[1 0 0 -1 0]`. However, `:A` will only be produced by the second pathway if a request for `:E` is made.
  - Make a copy of the material to specify an alternate pathway. This will require specifying parameters for the copied material throughout the network.
- Capacity limitations on shared feedstock inventory among producer nodes (e.g., shared inventory tanks) cannot be enforced directly. This is because the shared inventory is its own node and feeds the inventory holding area in the producer node. Thus the total inventory is the inventory at the inventory node plus the inventory positions at the producers. Capacity limitations must be enforced manually via the reorder actions. Potential fixes: (requires changing the package code)
  - Make the inventory capacity dynamic (a function of the producer inventory holding sites).
  - Split production site from feedstock inventory into two nodes. Requires updating model logic and behavior.
- If a `producer` can produce more than 1 material, it is possible to produce all materials it is capable of producing simultaneously. This does not account for resource constraints (e.g., single reactor can only do reaction 1 or reaction 2, but not both simultaneously). However, these can be enforced manually with the reorder actions. Potential fixes: (requires changing the package code)
  - Drop inventory capacities to 0 when the production equipment is occupied. Requires modeling each production unit as its own node.
  - Develop a production model (perhaps based on the Resource-Task Network paradigm)

## Inventory replenishment policies

At each iteration in the simulation, an `action` can be provided to the system, which consists of the replenishment orders placed on every link in the supply network. This `action` must be of type `Vector{Real}` and must be `nonnegative` of the form: `[Edge1_Material1, Edge1_Material2, ..., Edge1_MaterialM, Edge2_Material1, ..., Edge2_MaterialM, ..., EdgeE_Material1, ..., EdgeE_MaterialM]`, where the ordering in the edges is given by `edges(env.network)` and the ordering in the materials by `env.materials`.

An `action` vector can be visualized as a `DataFrame` using `show_action(action, env::SupplyChainEnv)`.

The function `reorder_policy` can be used to implement an inventory reorder policy at each node based its inventory position. Reorder quantities are placed to the node's priority supplier. The two most common policies used in industry are the `(s,S)` and `(r,Q)` [policies](https://smartcorp.com/inventory-control/inventory-control-policies-software/).

The `reorder_policy` takes the following inputs and returns an `action` vector.
- `env::SupplyChainEnv`: inventory management environment
- `reorder_point::Dict`: the `s` or `r` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_param::Dict`: the `S` or `Q` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_type::Union{Symbol, Dict}`: `:rQ` for an `(r,Q)` policy, or `:sS` for an `(s,S)` policy. If passing a `Dict`, the policy type should be specified for each node (keys).
- `review_period::Union{Int, StepRange, Vector, Dict}`: number of periods between each inventory review (Default = `1` for continuous review.). If a `StepRange` or `Vector` is used, the `review_period` indicates which periods the review is performed on. If a `Dict` is used, the review period should be specified for each `(node, material)` `Tuple` (keys). The values of this `Dict` can be either `Int`, `StepRange`, or `Vector`. Any missing `(node, material)` key will be assigned a default value of 1.
- `min_order_qty::Union{Real, Dict}`: minimum order quantity (MOQ) at each supply node. If a `Dict` is passed, the MOQ should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of 0. 

## Model Inputs

The supply network topology must be mapped on a network graph using [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl). The system parameters are stored in the network's metadata.

### Node-specific

`Producers` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: (*only when the node has at least 1 supplier*) `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.options[:reallocate] == true`).
- `:production_capacity::Dict`: maximum production capacity for each material (`keys`).

`Distributors` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.options[:reallocate] == true`).

`Markets` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.options[:reallocate] == true`).
- `:demand_distribution::Dict`: probability distributions from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) for the market demands for each material (`keys`). For deterministic demand, instead of using a probability distribution, use `D::Number`.
- `:demand_period::Dict`: mean number of periods between demand arrivals for each material (`keys`)
- `:demand_sequence::Dict`: a user specified `Vector` of market demand for each material (`keys`). When a nonzero `Vector` is provided, the `demand_distribution` and `demand_period` parameters are ignored.
- `:sales_price::Dict`: market sales price for each material (`keys`)
- `:stockout_penalty::Dict`: unit penalty for unsatisfied market demand for each material (`keys`)

### Edge-specific

All edges have the following fields in their metadata:
- `:sales_price::Dict`: unit sales price for inventory sent on that edge (from supplier to receiver) for each material (`keys`)
- `:transportation_cost::Dict`: unit transportation cost for shipped inventory for each material (`keys`)
- `:pipeline_holding_cost::Dict`: unit holding cost per period for inventory in-transit for each material (`keys`)
- `:lead_time::Distribution{Univariate, Discrete}`: probability distributions from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) for the lead times for each material (`keys`) on that edge. For deterministic lead times, instead of using a probability distribution, use `L::Number`.

### General Network

The graph metadata should have the following fields in its metadata:
- `:materials::Vector` with a list of all materials in the system.
- `:bill_of_materials::Matrix`: bill of materials indicating the production recipies for the materials in the system. The row numbers correspond to the input materials and the column numbers to the output materials. The numbering matches that of the `materials` vector. The magnitude of each element is proportional to the production of one unit of output material. Each element can have one of three types of values:
  - `zero`: input not involved in production of output.
  - `negative number`: input is consumed in the production of output.
  - `positive number`: input is a co-product of the output.

## Model Output

A `SupplyChainEnv` has the following fields:
- `network::MetaDiGraph`: Supply Chain Network (metagraph)
- `markets::Array`: list of market nodes
- `producers::Array`: list of producer nodes
- `materials::Array`: list of all material (material) names (strings)
- `inventory_on_hand::DataFrame`: timeseries with on hand inventories @ each node.
- `inventory_level::DataFrame`: timeseries with inventory level @ each node (on-hand minus backlog, if backlogging is allowed)
- `inventory_pipeline::DataFrame`: timeseries with pipeline inventories on each arc.
- `inventory_position::DataFrame`: timeseries with inventory positions @ each node (inventory level + placed replenishments).
- `replenishments::DataFrame`: timeseries Replenishment orders placed on each edge at the end of each period
- `shipments::DataFrame`: current shipments and time to arrival for each node
- `demand::DataFrame`: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market
- `profit::DataFrame`: timeseries with profit at each node
- `reward::Float64`: reward in the system (used for RL)
- `period::Int`: period in the simulation
- `num_periods::Int`: number of periods in the simulation
- `discount::Float64`: time discount factor (interest rate)
- `options::Dict`: simulation options:
  - `backlog::Bool`: backlogging allowed if `true`; otherwise, unfulfilled demand is lost sales
  - `reallocate::Bool`: the system try to reallocate requests if they cannot be satisfied if `true`; otherwise, no reallocation is attempted.
  - `evaluate_profit::Bool`: the simulation will calculate the proft at each node if `true` and save the results in `env.profit`.
- `seed::Int`: random seed

## Examples

### Example 1

This example is for a 100 period simulation of a supply network with one plant (node 1) that supplies a retailer (node 3), with stochastic demand for product `:A`. Node 3, has an alternate supplier, which is a distribution center (node 2). Node 3 prefers replenishing from the plant, which has a lower lead time. A `(s,S)` reorder policy is used at the retailer. When the on-hand level for material `:A` is depleted at the plant, the plant begins transforming raw material `:B` into `:A`. There is limited raw material supply at the plant. When the raw material stocks-out, node 3 switches to node 2 for its supply.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex1.jl).*

![](examples/figs/ex1_profit.png)
![](examples/figs/ex1_position.png)

### Example 2

This example is for a 100 period simulation of a supply network with one warehouse (node 1) that supplies a retailer (node 2) with stochastic demand for product `:A`. A `(r,Q)` reorder policy is used at the retailer every 25 periods.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex2.jl).*

![](examples/figs/ex2_inventory.png)

## Contact

**Author**: Hector D. Perez\
**Position**: Ph. D. Candidate @ Carnegie Mellon University\
**Email**: hdperez@cmu.edu\
**Year**: 2021
