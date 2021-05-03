# Sequence of Events

The following sequence of events occurs in each period of the simulation:
1. Start period.
2. Place inventory replenishment orders at each node. These are limited to the available production capacity (supplier is a producer) or available inventory. If `reallocate = true`, then any amount that cannot be satisfied is reallocated to the next priority supplier.
   - Distributors ship inventory.
   - Producers attempt to satisfy the requested amount in the following order: 1) using any on-hand inventory, 2) using any uncommitted scheduled production (only relevant when there is co-production), and 3) manufacturing material. Materials scheduled for production have a production lead time. Production costs are incurred at the start of production.
   - Producers send orders that have completed (after the production lead time).
4. Receive inventory that has arrived at each node (after the lead time has transpired).
5. Pay suppliers for inventory received.
6. Pay shipper for inventory shipped.
7. Market demand occurs after tossing a weighted coin with the probability of demand occurring defined by the `demand_frequency`.
8. Demand (including any backlog if `backlog = true`) is fulfilled up to available inventory at the markets.
9. Unfulfilled demand is penalized and backlogged (if `backlog = true`).
10. Each node pays a holding cost and a transportation cost for on-hand inventory and in-transit inventory at each period.
