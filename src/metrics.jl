"""
    calculate_service_measures!(env::SupplyChainEnv; window::Tuple=(0,Inf), metrics::Tuple=(:service_levels,:fill_rates))

Calculate service metrics specified in `metric` for each arc and material in the simulation.
    `window` sets the window of periods to apply in the service measure calculation.

NOTE:
- If reallocation occurs, any fulfilled reallocated order is not counted
- Node specific metrics are not pooled
"""
function calculate_service_measures!(env::SupplyChainEnv; window::Tuple=(0,env.period), metrics::Tuple=(:service_levels,:fill_rates))
    @assert window[1] <= window[2] "`window` is not valid (lb > ub)."
    if window[1] >= env.period
        @warn "Service measure window outside of simulation time. `calculate_service_measures! not called."
        return
    end
    # @assert window[1] <= env.period "Service measure window outside of simulation."
    if window == (0,Inf)
        orders_df = env.orders
        fulfillments_df = env.fulfillments
    else
        orders_df = filter([:created,:due] => (c,d) -> window[1] <= c && d <= window[2], env.orders, view = true)
        id_set = Set(orders_df.id)
        fulfillments_df = filter(:id => in(id_set), env.fulfillments, view=true)
    end
    if :fill_rates in metrics
        env.metrics[:fill_rates] = calculate_fill_rates(orders_df, fulfillments_df)
    end
    if :service_levels in metrics
        env.metrics[:service_levels] = calculate_service_levels(orders_df)
    end
end

function calculate_fill_rates(orders_df::AbstractDataFrame, fulfillments_df::AbstractDataFrame)
    #filter out late orders and lost_sales
    due_dict = Dict(orders_df.id .=> orders_df.due)
    @chain fulfillments_df begin
        subset( #Only keep fulfillments specified by fulfillment type and remove late fulfillments
            :type => ByRow(==(:sent)),
            [:period, :id] => ByRow((t,id) -> t <= due_dict[id]),
            view = true
        )
        groupby([:id,:material,:arc])
        combine(:amount => sum => :filled) #sum deliveries on each order
        leftjoin(orders_df, _, on = [:id,:material,:arc]) #NOTE: will exclude reallocated orders (have a different arc on the fulfillments df)
        transform([:filled, :amount] => ByRow((f,a) -> f/a) => :fill) #percentage filled (fill)
        coalesce.(0) #replace missings with 0 (not filled)
        groupby([:material,:arc])
        combine(:fill => mean => :fill_rate)
    end
end

calculate_service_levels(orders_df::AbstractDataFrame) = 
    @chain orders_df begin
        transform(
            [:fulfilled, :due] => ByRow((f,d) -> ismissing(f) || f == :lost_sale ? false : f <= d) => :on_time
        )
        groupby([:material,:arc])
        combine(
            :on_time => mean => :service_level
        )
    end