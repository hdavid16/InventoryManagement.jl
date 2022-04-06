# Inventory replenishment policies

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
