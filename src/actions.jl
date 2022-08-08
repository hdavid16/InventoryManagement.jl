"""
    (x::SupplyChainEnv, action::Vector{T} where T <: Real)

Convert a replenishment order vector into a NamedArray indicating how much of
each material (rows) is being requested on each arc (columns).
"""
function show_action(x::SupplyChainEnv, action::Vector{T} where T <: Real)
    mats = x.materials
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    action_named_array = NamedArray(
        reshape(action, (length(mats), length(arcs))),
        (x.materials, arcs),
        (:material, :arc)
    )
    
    return action_named_array
end

"""
    (x::SupplyChainEnv)(action::Vector{T} where T <: Real)

Apply an `action` (replenishment requests) on the `SupplyChainEnv` and step
forward one simulation period.
"""
function (x::SupplyChainEnv)(action::Vector{T} where T <: Real)
    #validate action input and reshape to matrix
    @assert all(action .>= 0) "Reorder actions cannot be negative."
    @assert length(action) == length(x.materials)*ne(x.network) "Reorder action vector must have length num_products * num_edges."
    act = show_action(x, action)
    #increase period counter
    x.period += 1
    #intialize next period on-hand and pipeline inventories with previous inventories
    initialize_inventories!(x)
    #move active shipments forward one period
    x.shipments.lead .-= 1
    #decrease time until due date for active orders by one period
    x.open_orders.due .-= 1
    #receive incoming orders
    update_shipments!(x)
    #place requests
    place_orders!(x, act)
    #discard any excess inventory
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)
    #markets open and demand occurs
    simulate_markets!(x)
    #update inventories at each node and echelon
    update_inventories!(x)
    #calculate profit at each node
    if x.options[:evaluate_profit]
        calculate_profit!(x)
        x.reward = @chain x.profit begin
            @rsubset(:period == x.period)
            sum(_.value; init=0)
        end
        x.reward = sum(@rsubset(x.profit, :period == x.period, view=true).value) #update reward (current profit). NOTE: is this ok for RL?
    end
end
