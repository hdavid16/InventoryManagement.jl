# Model

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

## Model Inputs

The supply network topology must be mapped on a network graph using [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl). The system parameters are stored in the network's metadata.

### Node-specific

`Producers` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: (*only when the node has at least 1 supplier*) `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).
- `:production_cost::Dict`: unit production cost for each material (`keys`)
- `:production_capacity::Dict`: maximum production capacity for each material (`keys`).
- `:production_time::Dict`: production lead time for each material (`keys`).

`Distributors` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).

`Markets` will have the following fields in their node metadata:
- `:initial_inventory::Dict`: initial inventory for each material (`keys`)
- `:inventory_capacity::Dict`: maximum inventory for each material (`keys`)
- `:holding_cost::Dict`: unit holding cost for each material (`keys`)
- `:supplier_priority::Dict`: `Vector` of supplier priorities (from high to low) for each material (`keys`). When a request cannot be fulfilled due to insufficient productio capacity or on-hand inventory, the system will try to reallocate it to the supplier that is next in line on the priority list (if `env.reallocate == true`).
- `:demand_distribution::Dict`: probability distributions from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) for the market demands for each material (`keys`). For deterministic demand, instead of using a probability distribution, use `[D]` where `D` is a `Number`.
- `:demand_frequency::Dict`: probability that demand will occur (value between `0.0` and `1.0`) for each material (`keys`)
- `:demand_sequence::Dict`: a user specified `Vector` of market demand for each material (`keys`). When a nonzero `Vector` is provided, the `demand_distribution` and `demand_frequency` parameters are ignored.
- `:sales_price::Dict`: market sales price for each material (`keys`)
- `:demand_penalty::Dict`: unit penalty for unsatisfied market demand for each material (`keys`)

### Edge-specific

All edges have the following fields in their metadata:
- `:sales_price::Dict`: unit sales price for inventory sent on that edge (from supplier to receiver) for each material (`keys`)
- `:transportation_cost::Dict`: unit transportation cost per period for inventory in-transit for each material (`keys`)
- `:lead_time::Distribution{Univariate, Discrete}`: probability distributions from [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) for the lead times for each material (`keys`) on that edge. For deterministic lead times, instead of using a probability distribution, use `[L]` where `L` is a `Number`.

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
- `distributors::Array`: list of distribution nodes (excludes end distributors where markets exist)
- `materials::Array`: list of all material (material) names (strings)
- `bill_of_materials::Matrix` square matrix with BOM (rows = input materials, cols = output materials; indices follow materials list; positive value is a co-product, negative is a feedstock)
- `inv_on_hand::DataFrame`: timeseries On Hand Inventory @ each node at the end of each period
- `inv_pipeline::DataFramet`: timeseries Pipeline Inventory on each edge at the end of each period
- `inv_position::DataFrame`: timeseries Inventory Position for each node at the end of each period
- `replenishments::DataFrame`: timeseries Replenishment orders placed on each edge at the end of each period
- `shipments::DataFrame`: current shipments and time to arrival for each node
- `production::DataFrame`: current material production committed to an edge and lead time to ship. Note: byproducts are scheduled to go to the producing node `n` (edge `(n,n)`).
- `demand::DataFrame`: timeseries with realization of demand, sold units, unfulfilled demand, and backlog at each market
- `profit::DataFrame`: timeseries with profit at each node
- `reward::Float64`: reward in the system (used for RL)
- `period::Int`: period in the simulation
- `num_periods::Int`: number of periods in the simulation
- `discount::Float64`: time discount factor (interest rate)
- `backlog::Bool`: backlogging allowed if `true`; otherwise, unfulfilled demand is lost sales
- `reallocate::Bool`: the system try to reallocate requests if they cannot be satisfied if `true`; otherwise, no reallocation is attempted.
- `seed::Int`: random seed
