"""
    service_measures(env::SupplyChainEnv)

Calculate mean service level and fill rate for each node and each material in the simulation.

If there is more than 1 supplier and reallocation of requests occurs, the metric will be associated with the original supplier.
"""
function service_measures(env::SupplyChainEnv; review_period::Union{Int, StepRange, Vector, Dict} = 1)
    #filter out times with no demand
    demand = filter(i -> i.quantity > 0, env.demand) 
    #generate a column with the supplier for each order
    transform!(demand,
        :arc => ByRow(first) => :supplier
    )
    #aggregate downstream demand and determine for each period, material, and supplier, the percentage of the order that is filled and whether the order is fully filled
    pooled_demand = combine(
        groupby(demand, [:period, :material, :supplier]),
        [:quantity, :fulfilled] => 
            ((q,f) -> (sum(q), sum(f)) |> 
                sqf -> (fraction_filled = sqf[2]/sqf[1], order_filled = sqf[1] == sqf[2]))
            => AsTable
    )
    #get average of fraction_filled and order_filled for each material and supplier to get FR and CSL respectively
    service_measures = combine(
        groupby(pooled_demand, [:material, :supplier]),
        :fraction_filled => mean => :fill_rate,
        :order_filled => mean => :service_level
    )
    
    return service_measures
end