"""
    inventory_balance(env::SupplyChainEnv, node::Int, material::Material, inventory_type::Symbol = :on_hand)

Create a `DataFrame` logging the following in each period:
1) opening inventory (for the specified inventory type; options: `:on_hand`, `:level`, or `:position`)
2) new external demand
3) fulfillment of external demand
4) new internal demand
5) fulfillment of internal demand
6) arrivals from upstream suppliers 
7) unfulfilled demand
"""
function inventory_balance(env::SupplyChainEnv, node::Int, material::Material, inventory_type::Symbol = :on_hand)
    #on_hand inventory
    balance = opening_inventory(env,node,material,inventory_type)
    #external demand
    if node in env.markets
        external_demand!(balance,env,node,material)
        external_fulfillment!(balance,env,node,material)
    else
        insertcols!(balance, :ext_demand => 0)
        insertcols!(balance, :ext_backlog => 0)
    end
    #internal demand
    successors = outneighbors(env.network,node)
    if !isempty(successors)
        internal_demand!(balance,env,node,material)
        internal_fulfillment!(balance,env,node,material)
    else
        insertcols!(balance, :int_demand => 0)
        insertcols!(balance, :int_backlog => 0)
    end
    #replenishment arrivals
    predecessors = inneighbors(env.network,node)
    if !isempty(predecessors)
        replenishments!(balance,env,node,material)
    else
        insertcols!(balance, :arrivals => 0)
    end
    #unfulfilled demand
    unfulfilled!(balance,env,node,material)
end

"""
    opening_inventory(env::SupplyChainEnv, node::Int, material::Material, inventory_type::Symbol)

Create a `DataFrame` with the inventory time series.
"""
function opening_inventory(env::SupplyChainEnv, node::Int, material::Material, inventory_type::Symbol)
    @assert inventory_type in [:on_hand, :level, :position] "Invalid `inventory_type` passed. Valid options are `:on_hand`, `:level`, or `:position`"
    @chain env.inventory begin
        subset( #filter inventory type, node, material
            :type => ByRow(==(inventory_type)), 
            :location => ByRow(==(node)),
            :material => ByRow(==(material)),
            view = true
        )
        select(:period, :material, :amount => :starting)
        transform!(:period => ByRow(p -> p + 1) => :period) #shift time
    end
end

"""
    external_demand!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material)

Append the external demand at each period to the `balance` DataFrame.
"""
external_demand!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) =
    @chain env.orders begin
        subset( #get market orders for that material at that node
            :arc => ByRow(==((node,:market))),
            :material => ByRow(==(material)),
            view = true
        )
        groupby(:created)
        combine(:amount => sum => :ext_demand) #add all orders on same period
        leftjoin!(balance, _, on = :period => :created) #merge with balance
    end

"""
    external_fulfillment!(balance::DataFrame,env::SupplyChainEnv, node::Int, material::Material)

Append the external demand fulfillment at each period to the `balance` DataFrame.
"""
external_fulfillment!(balance::DataFrame,env::SupplyChainEnv, node::Int, material::Material) = 
    @chain env.fulfillments begin
        subset( #get order fulfillments
            :type => ByRow(==(:sent)), #fulfilled when sent
            :arc => ByRow(==((node,:market))),
            :material => ByRow(==(material)),
            view = true
        )
        groupby(:period)
        combine(:amount => sum => :ext_fulfill)
        leftjoin!(balance, _, on = :period)
    end

"""
    internal_demand!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material)

Append the internal demand at each period to the `balance` DataFrame.
"""
internal_demand!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) =
    @chain env.orders begin
        subset( #get downstream orders
            :arc => ByRow(a -> a[1] == node && a[2] != :market),
            :material => ByRow(==(material)),
            view = true
        )
        groupby(:created)
        combine(:amount => sum => :int_demand) #add all orders on same period
        leftjoin!(balance, _, on = :period => :created) #merge with balance
    end

"""
    internal_fulfillment!(balance::DataFrame,env::SupplyChainEnv, node::Int, material::Material)

Append the internal demand fulfillment at each period to the `balance` DataFrame.
"""
internal_fulfillment!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) = 
    @chain env.fulfillments begin
        subset( #get order fulfillments
            :type => ByRow(==(:sent)), #fulfilled when sent
            :arc => ByRow(a -> a[1] == node && a[2] != :market),
            :material => ByRow(==(material)),
            view = true
        )
        groupby(:period)
        combine(:amount => sum => :int_fulfill)
        leftjoin!(balance, _, on = :period)
    end

"""
    replenishments!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material)

Append the arrived replenishments at each period to the `balance` DataFrame
"""
replenishments!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) = 
    @chain env.fulfillments begin
        subset( #get arrivals
            :type => ByRow(==(:delivered)),
            :arc => ByRow(a -> a[2] == node),
            :material => ByRow(==(material)),
            view = true
        )
        groupby(:period)
        combine(:amount => sum => :arrivals) #add all orders on same period
        leftjoin!(balance, _, on = :period) #merge with balance
    end

"""
    unfulfilled!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) = 

Append total unfulfilled demand at each period to the `balance` DataFrame
"""
unfulfilled!(balance::DataFrame, env::SupplyChainEnv, node::Int, material::Material) = 
    @chain env.inventory begin
        subset( #filter inventory type, node, material
            :type => ByRow(==(:unfulfilled)), 
            :location => ByRow(==(node)),
            :material => ByRow(==(material)),
            view = true
        )
        select(:period, :material, :amount => :unfulfilled)
        leftjoin!(balance, _, on = [:period, :material])
    end