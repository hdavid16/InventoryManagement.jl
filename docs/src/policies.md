# Inventory replenishment policies

At each iteration in the simulation, an `action` can be provided to the system, which consists of the replenishment orders placed on every link in the supply network. This `action` must be of type `Vector{Real}` and must be `nonnegative` of the form: `[Arc_1_Material_1, Arc_1_Material_2, ..., Arc_1Material_M, Arc_2_Material_1, ..., Arc_2_Material_M, ..., Arc_A_Material_1, ..., Arc_A_Material_M]`, where the ordering in the arcs is given by `edges(SupplyChainEnv.network)` and the ordering in the materials by `SupplyChainEnv.materials`.

An `action` vector can be visualized as a `NamedArray` using `show_action(SupplyChainEnv, action)`:

```julia
material ╲ arc │ :Arc_1  :Arc_2 ... :Arc_A
───────────────┼──────────────────────────
:Material_1    │  
:Material_2    │  
...            │
:Material_M    │  
```

The function `reorder_policy` can be used to implement an inventory reorder policy at each node based its inventory position or echelon stock. Reorder quantities are placed to the node's priority supplier. The reorder policy is applied for each `material` at each `node` in reverse [topological order](https://en.wikipedia.org/wiki/Topological_sorting). This allows upstream nodes to determine their reorder quantities with information about the reorder quantities placed by their successors (relevant for `producer` nodes to ensure that raw material replenishments are synced with production orders). The two most common policies used in industry are the `(s,S)` and `(r,Q)` [policies](https://smartcorp.com/inventory-control/inventory-control-policies-software/).

The `reorder_policy` takes the following inputs and returns an `action` vector.
- `env::SupplyChainEnv`: inventory management environment
- `reorder_point::Dict`: the `s` or `r` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_param::Dict`: the `S` or `Q` parameter in each node for each material in the system. The `keys` are of the form `(node, material)`.
- `policy_type::Union{Symbol, Dict}`: `:rQ` for an `(r,Q)` policy, or `:sS` for an `(s,S)` policy. If passing a `Dict`, the policy type should be specified for each node (`keys`).
- `review_period::Union{Int, AbstractRange, Vector, Dict}`: number of periods between each inventory review (Default = `1` for continuous review.). If a `AbstractRange` or `Vector` is used, the `review_period` indicates which periods the review is performed on. If a `Dict` is used, the review period should be specified for each `(node, material)` `Tuple` (`keys`). The values of this `Dict` can be either `Int`, `AbstractRange`, or `Vector`. Any missing `(node, material)` key will be assigned a default value of 1.
- `min_order_qty::Union{Real, Dict}`: minimum order quantity (MOQ) at each supply node. If a `Dict` is passed, the MOQ should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of 0. 
- `order_multiples::Union{Real, Dict}`: size increments for each order (default is -1, which means no constraint on order sizes). If a `Dict` is passed, the order multiples should be specified for each `(node, material)` `Tuple` (keys). The values should be `Real`. Any missing key will be assigned a default value of -1 (meaning no order multiples enforced).
- `adjust_expected_consumption::Bool`: indicator if the reorder point should be increased (temporarilly) at a `producer` node by the expected raw material consumption for an expected incoming production order.