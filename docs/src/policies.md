# Inventory replenishment policies

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
