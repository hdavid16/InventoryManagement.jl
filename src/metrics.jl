"""
    service_measures(env::SupplyChainEnv)

Calculate mean service level and fill rate for each node and each material in the simulation.

NOTE: Not suported when there is more than 1 supplier and reallocation is allowed.
      CSL is currently on a per-period basis, not a cycle basis
"""
function service_measures(env::SupplyChainEnv; review_period::Union{Int, StepRange, Vector, Dict} = 1)
    #merge demand and replenishment tables
    demand_filt = filter(i -> i.quantity > 0, env.orders) #filter out times with no demand at markets
    replenish_filt = filter(i -> i.quantity > 0, env.orders)
    select!(replenish_filt, #convert replenishment table into same format as demand table for merging
        :period,
        :arc => ByRow(first) => :node,
        :material,
        :quantity => :demand,
        :fulfilled => :sold,
        :unfulfilled
    )
    orders = vcat(demand_filt, replenish_filt) #merge two tables

    #convert orders of products at producer nodes into orders of raw materials at those nodes
    # bom = env.bill_of_materials
    for row in eachrow(filter(i -> i.node in env.producers, orders))
        bom = get_prop(env.network, row.node, :bill_of_materials)
        col = row.material
        row_mats = names(bom,1)
        for id in findall(i -> i < 0, bom[:,col])
            push!(orders, (row.period, row.node, row_mats[id], -row.quantity*bom[id,col], -row.sold*bom[id,col], -row.unfulfilled*bom[id,col]))
        end
    end

    #convert review period into proper format (dictionary by node)
    if review_period isa Int
        review_period = 1:review_period:env.num_periods
    end
    if !isa(review_period, Dict)
        review_period = Dict((n,mat) => review_period for n in vertices(env.network), mat in env.materials)
    end
    cycle = Dict(
        key => 
            Dict(
                i => 
                    Interval{:closed, :open}(review_period[key][i],review_period[key][i+1])
                for i in 1:length(review_period[key])-1
            ) 
        for key in keys(review_period)
    )

    #calculate CSL
    orders.cycle = Any[nothing for _ in eachrow(orders)] #initialize column to store what cycle a row belongs to
    orders_filt = filter(
        i -> !isnothing(i.cycle),
        combine(
            groupby(orders, [:node, :material]),
            [:period,:node,:material] => ByRow((t,n,mat) -> findfirst(x -> t in x, cycle[n,mat])) => :cycle,
            :unfulfilled
        )
    )
    stockouts = combine(
        groupby(orders_filt, [:cycle, :node, :material]), 
        :unfulfilled => !iszero ∘ sum => :stockout #count the if a stockout occurs in each period
    )
    CSL_stack = combine(
        groupby(stockouts, [:node, :material]), 
        :stockout => (x -> 1 - mean(x)) => :CycleServiceLevel
    )
    CSL = sort(
        unstack(CSL_stack, :node, :material, :CycleServiceLevel), 
        :node
    )

    #calculate FR
    FR_stack = combine(
        groupby(orders, [:node, :material]), 
        [:demand, :sold] => ((x,y) -> sum(y)/sum(x)) => :FillRate
    )
    FR = sort(
        unstack(FR_stack, :node, :material, :FillRate), 
        :node
    )
    
    return FR, CSL
end