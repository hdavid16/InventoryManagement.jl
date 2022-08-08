# Sequence of Events

The sequence of events in each period of the simulation is patterned after that of the [News Vendor Problem](https://optimization.cbe.cornell.edu/index.php?title=Newsvendor_problem):
1. Start period.
2. Receive any incoming shipments from upstream nodes.
3. Place inventory replenishment orders at each node by traversing the supply network downstream (using [topological sorting](https://en.wikipedia.org/wiki/Topological_sorting)).
  - If `backlog = true`, the previous period's backlog is added to the replenishment order. 
  - Each supplier attempts to fulfill downstream orders via its on-hand inventory. `Producer` nodes fulfill production requests via material production (if there is sufficient `production capacity` and `raw material inventory`). 
  - If `reallocate = true`, then any amount that cannot be satisfied is reallocated to the next supplier in the `supplier priority` list (the lowest priority supplier will reallocate back to the highest priority supplier). 
  - Accepted replenishment orders are immediately shipped with a lead time sampled from the specified distribution. For `distributor` nodes, the lead time is the in-transit (transportation) time between `distributor` nodes. For `producer` nodes, the lead time is the plant production time. 
  - If the lead time is 0, the stock will be immediately available to the requesting node so that it can be used to fulfill downstream orders as they arrive (possible due to the topological sorting).
4. Market (external) demand for each material occurs after tossing a weighted coin with the probability of demand occurring defined by the `demand_frequency` (likelihood of demand occuring in each period; inverse of the average number of periods between positive external demands).
  - For example, if `demand_requency = 0.5`, there is a `50%` chance of having positive demand, or once every 2 days on average.
5. Demand (including any backlog if `backlog = true`) is fulfilled up to available inventory at the `market` nodes. Make-to-order nodes trigger production requests.
6. Unfulfilled demand is backlogged (if `backlog = true`).
7. Accounts for each node are generated:
  - Accounts payable: invoice for fulfilled replenishment orders (payable to suppliers), invoice for delivered replenishment orders (payable to third-party shipper), pipeline inventory holding cost for in-transit inventory (cost to requestor), on-hand inventory holding cost, penalties for unfulfilled demand (cost to supplier).
  - Accounts receivable: sales for internal and external demand.
