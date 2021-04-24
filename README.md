# InventoryManagement.jl:

*Discrete-time simulation environment for Inventory Management in Supply Networks.*

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
9. [Model Output](#model-output)
10. [Examples](#examples)
    - [Example #1: alternate suppliers and continuous review (s,S) policy](#example-1)
    - [Example #2: unlimited upstream supply and periodic (r,Q) policy](#example-2)
    - [Example #3: make-to-stock plant with market demand](#example-3)
    - [Example #4: chemical production system with co-production and material recycle](#example-4)
11. [Contact](#contact)

## Overview

*InventoryManagement.jl* allows modeling a multi-period multi-product supply network. A supply network can be constructed using the following node types:
- `Producers`: Nodes where inventory transformation takes place (e.g., intermediates or final materials are produced). Reactive systems, including those with co-products, can be modelled using [Bills of Materials](https://en.wikipedia.org/wiki/Bill_of_materials) (see [Model Inputs section](#graph-specific)).
- `Distributors`: Intermediate nodes where inventory is stored and distributed (e.g., distribution centers).
- `Markets`: Nodes where end-customers place final product orders (i.e., retailer). These are the last (sink) nodes in the network.

The simplest network that can be modeled is one with a single market with one producer or distributor. However, more complex systems can be modelled as well.

When defining a supply network, a `SupplyChainEnv` object is created based on system inputs and network structure. This object can then be used to execute a simulation of the inventory dynamics. During a simulation, stochastic demand at each of the markets can occur for each of the materials in each period. When product demand occurs at the market, sales are made based on available inventory. Any unfulfilled demand is either backlogged or considered a lost sale depending on the system definition. If no action is taken duirng the simulation, the inventory levels will eventually be depleted. To avoid this from happening, a decision-maker can interact with the system in each period by making inventory replenishment decisions (refered to as `actions`). Lead times for in-transit inventory as well as production lead times are accounted for in the simulation. Transportation lead times can be modelled stochastically to account for lead time uncertainty. From a service time perspective, demand at market nodes has zero service time, whereas non-market nodes have service time equal to the production lead time + transportation lead time.

The `SupplyChainEnv` can also potentially be used in conjunction with [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl) to train a Reinforcement Learning `agent`.

This package generalizes and extends and the inventory management environment available in [OR-Gym](https://github.com/hubbs5/or-gym).

## Dependencies

*InventoryManagement.jl* mainly relies on the following packages:
- [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl): Define supply network structure and specify node- and edge-specific parameters.
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl): Tabulate results.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Define probability distributions for the lead times in between nodes and the demands at the market nodes.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```julia
pkg> add https://github.com/hdavid16/InventoryManagement.jl
```

## Sequence of Events

The following sequence of events occurs in each period of the simulation:
1. Start period.
2. Place inventory replenishment orders at each node. These are limited to the available production capacity (supplier is a producer) or available inventory. If `reallocate = true`, then any amount that cannot be satisfied is reallocated to the next priority supplier.
   - Distributors ship inventory.
   - Producers attempt to satisfy order with any on-hand inventory and then manufacture materials for any amount remaining (a production lead time begins). Production costs are incurred at the start of production.
   - Producers send orders that have completed (after the production lead time).
4. Receive inventory that has arrived at each node (after the lead time has transpired).
5. Pay suppliers for inventory received.
6. Pay shipper for inventory shipped.
7. Market demand occurs after tossing a weighted coin with the probability of demand occurring defined by the `demand_frequency`.
8. Demand (including any backlog if `backlog = true`) is fulfilled up to available inventory at the markets.
9. Unfulfilled demand is penalized and backlogged (if `backlog = true`).
10. Each node pays a holding cost and a transportation cost for on-hand inventory and in-transit inventory at each period.

## Model Assumptions

The following assumptions hold in the current implementation, but can be modified in future releases.

- `Producers` produce material on demand ([make-to-order](https://en.wikipedia.org/wiki/Build_to_order) policy).
- `Producers` can hold inventory. Downstream replenishment orders are fulfilled first with any on-hand inventory, and then via production only after there is no on-hand inventory left. 
- Replenishment orders can only be satisfied with current on-hand inventory or available production capacity.
- Commited production orders count towards the inventory position of the downstream node.
- Backlogging is only allowed at the `Markets`, it is not allowed for inventory replenishment decisions.
- Transportation costs are paid to a third party (not a node in the network).

## Model Limitations

The following features are not currently supported:

- Inventory capacity limits at the nodes is not currently supported.
- `Producers` do not operate under a [make-to-stock](https://en.wikipedia.org/wiki/Build_to_stock) policy since any material produced gets shipped downstream. However, this can be accomodated by adding a dumby node downstream of the `producer` that holds inventory produced by the `producer` with zero lead time in between the nodes. Thus, using a proper reorder policy, the `producer` can act as a `make-to-stock` system that pushes inventory to the inventory holding node.
- `Producers` do not have market demand. However, this can be modelled by adding a `market` node downstream of the `producer` with zero lead time in between the nodes.
- Alternate bills of materials (see [Model Inputs](#graph-specific)) for the same material are not currently supported. This is particularly relevant for chemical systems. However, the following workarounds can be done:
  - If the alternate reaction pathway has a byproduct, then the main product can be included as a co-product in the bill of materials of the byproduct. For example: A system with 5 materials (`:A - :E`) can have two ways to produce `:A`, `:B + :C -> :A` and `:D -> :A + :E`. The column for material `:A` can have the bill of material: `[0 -1 -1 0 0]`. The column for material `:E` can have the bill of materials: `[1 0 0 -1 0]`. However, `:A` will only be produced by the second pathway if a request for `:E` is made.
  - Make a copy of the material to specify an alternate pathway. This will require specifying parameters for the copied material throughout the network.

## Inventory replenishment policies

At each iteration in the simulation, an `action` can be provided to the system, which consists of the replenishment orders placed on every link in the supply network. This `action` must be of type `Vector{Real}` and must be `nonnegative` of the form: `[Edge1_Material1, Edge1_Material2, ..., Edge1_MaterialM, Edge2_Material1, ..., Edge2_MaterialM, ..., EdgeE_Material1, ..., EdgeE_MaterialM]`, where the ordering in the edges is given by `edges(env.network)` and the ordering in the materials by `env.materials`.

An `action` vector can be visualized as a `DataFrame` using `show_action(action, env::SupplyChainEnv)`.

The function `reorder_policy` can be used to implement an inventory reorder policy. The two most common policies used in industry are the `(s,S)` and `(r,Q)` [policies](https://smartcorp.com/inventory-control/inventory-control-policies-software/).

The `reorder_policy` takes the following inputs and returns an `action` vector.
- `env::SupplyChainEnv`: inventory management environment
- `param1::Dict`: the `s` or `r` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `param2::Dict`: the `S` or `Q` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `level::Symbol`: `:position` if the policy is based on the node's inventory position, or `:on_hand` if the policy is based on the node's on-hand inventory level.
- `kind::Symbol`: `:rQ` for an `(r,Q)` policy, or `:sS` for an `(s,S)` policy
- `supplier_selection::Symbol`: evenly distribute reorder quantities among all suppliers if `:random`; otherwise (if `:priority`), assign reorder quantities based on supplier priority (e.g., if supplier 1 does not have enough capacity or inventory, then request as much as possible and then request any remaining amount from the next supplier, and so forth).
- `review_period::Int`: number of periods between each inventory review (Default = `1` for continuous review. *Note*: only relevant when `level = :on_hand`)

## Model Inputs

### Node-specific

`Producers` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: (*only when the node has at least 1 supplier*) `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).
- `:production_cost::Dict`: unit production cost for each material (`keys`)
- `:production_capacity::Dict`: maximum production capacity for each material (`keys`).
- `:production_time::Dict`: production lead time for each material (`keys`).

`Distributors` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).

`Markets` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).
- `:demand_distribution::Dict`: probability distributions from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) for the market demands for each material (`keys`). For deterministic demand, instead of using a probability distribution, use `[D]` where `D` is a `Number`.
- `:demand_frequency::Dict`: probability that demand will occur (value between `0.0` and `1.0`) for each material (`keys`)
- `:sales_price::Dict`: market sales price for each material (`keys`)
- `:demand_penalty::Dict`: unit penalty for unsatisfied market demand for each material (`keys`)

### Edge-specific

All edges have the following fields in their metadata:
- `:sales_price::Dict`: unit sales price for inventory sent on that edge (from supplier to receiver) for each material (`keys`)
- `:transportation_cost::Dict`: unit transportation cost per period for inventory in-transit for each material (`keys`)
- `:lead_time::Distribution{Univariate, Discrete}`: the lead time on each edge

### Graph-specific

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
- `distributors::Array`: list of distribution nodes (excludes end distributors where markets exist)
- `materials::Array`: list of all material (material) names (strings)
- `bill_of_materials::Matrix` square matrix with BOM (rows = input materials, cols = output materials; indices follow materials list; positive value is a co-product, negative is a feedstock)
- `inv_on_hand::DataFrame`: timeseries On Hand Inventory @ each node at the end of each period
- `inv_pipeline::DataFramet`: timeseries Pipeline Inventory on each edge at the end of each period
- `inv_position::DataFrame`: timeseries Inventory Position for each node at the end of each period
- `replenishments::DataFrame`: timeseries Replenishment orders placed on each edge at the end of each period
- `shipments::DataFrame`: current shipments and time to arrival for each node
- `production::DataFrame`: current material production committed to an edge and lead time to ship
- `demand::DataFrame`: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market
- `profit::DataFrame`: timeseries with profit at each node
- `reward::Float64`: reward in the system (used for RL)
- `period::Int`: period in the simulation
- `num_periods::Int`: number of periods in the simulation
- `discount::Float64`: time discount factor (interest rate)
- `backlog::Bool`: backlogging allowed if `true`; otherwise, unfulfilled demand is lost sales
- `reallocate::Bool`: the system try to reallocate requests if they cannot be satisfied if `true`; otherwise, no reallocation is attempted.
- `seed::Int`: random seed

## Examples

### Example 1

This example is for a 100 period simulation of a supply network with one plant (node 1) that supplies a retailer (node 3), with stochastic demand for product `:A`. Node 3, has an alternate supplier, which is a distribution center (node 2). Node 3 prefers replenishing from the plant, which has a lower lead time. A `(s,S)` reorder policy is used at the retailer. When the on-hand level for material `:A` is depleted at the plant, the plant begins transforming raw material `:B` into `:A`. There is limited raw material supply at the plant. When the raw material stocks-out, node 3 switches to node 2 for its supply.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex1.jl).*

![](examples/figs/ex1_profit.png)
![](examples/figs/ex1_position.png)

### Example 2

This example is for a 100 period simulation of a supply network with one warehouse (node 1) that supplies a retailer (node 2), with stochastic demand for product `:A`. A `(r,Q)` reorder policy is used at the retailer every 25 periods.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex2.jl).*

![](examples/figs/ex2_inventory.png)

### Example 3

### Example 4

## Contact

**Author**: Hector D. Perez\
**Position**: Ph. D. Candidate @ Carnegie Mellon University\
**Email**: hdperez@cmu.edu\
**Year**: 2021
