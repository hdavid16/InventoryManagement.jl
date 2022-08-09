"""
    calculate_service_measures!(env::SupplyChainEnv; window::Tuple=(0,Inf), fulfillment_type::Symbol = :delivered)
    calculate_service_measures(orders_df::DataFrame, fulfillments_df::DataFrame, fulfillment_type::Symbol = :delivered)

Calculate mean service level and fill rate along each arc and for each material in the simulation.
    `window` sets the window of periods to apply in the service measure calculation.

NOTE:
- If reallocation occurs, any fulfilled reallocated order is not counted
- Node specific metrics are not pooled
"""
function calculate_service_measures(orders_df::AbstractDataFrame, fulfillments_df::AbstractDataFrame, fulfillment_type::Symbol = :delivered)
    @assert fulfillment_type in [:delivered, :sent] "`fulfillment_type` must be either `:sent` (shipped) or `:delivered`."
    #filter out lost_sales
    due_dict = Dict(orders_df.id .=> orders_df.due)
    df = @chain fulfillments_df begin
        subset( #Only keep deliveires and remove late fulfillments
            :type => ByRow(==(fulfillment_type)),
            [:period, :id] => ByRow((t,id) -> t <= due_dict[id]),
            view = true
        )
        groupby([:id,:material,:arc])
        combine(:amount => sum => :fulfilled) #sum deliveries on each order
        leftjoin(orders_df, _, on = [:id,:material,:arc]) #NOTE: will exclude reallocated orders (have a different arc on the fulfillments df)
        transform(
            [:fulfilled, :amount] => ByRow(
                (f,a) -> (fill = f/a, service = f ≈ a)) => AsTable #percentage filled (fill); fulfuilled? (service)
        )
    end
    id_set = Set(subset(df, :service => ByRow(identity), view=true, skipmissing=true).id)
    @chain df begin
        filter( #filter out :consumption orders that were fulfilled from stock (e.g., MTO orders)
            [:id, :arc, :service] => (id,a,s) -> !(
                ismissing(s) &&
                a[2] == :consumption &&
                id in id_set
            ),
            _,
            view = true
        )
        coalesce.(0) #replace missings with 0 (not filled)
        groupby([:material,:arc])
        combine(
            :service => mean => :avg_service,
            :fill => mean => :avg_fill
        )
    end
end
function calculate_service_measures!(env::SupplyChainEnv; window::Tuple=(0,Inf), fulfillment_type::Symbol = :delivered)
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
        orders_df = filter(:created => t -> window[1] <= t <= window[2], env.orders, view = true)
        id_set = Set(orders_df.id)
        fulfillments_df = filter(:id => in(id_set), env.fulfillments, view=true)
    end
    env.metrics = calculate_service_measures(orders_df, fulfillments_df, fulfillment_type)
end