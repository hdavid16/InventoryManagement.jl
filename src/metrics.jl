"""
    calculate_service_measures!(env::SupplyChainEnv)
    calculate_service_measures(orders_df::DataFrame, fulfillments_df::DataFrame)

Calculate mean service level and fill rate along each arc and for each material in the simulation.

NOTE:
- If reallocation occurs, any fulfilled reallocated order is not counted
- Node specific metrics are not pooled
"""
function calculate_service_measures(orders_df::AbstractDataFrame, fulfillments_df::AbstractDataFrame)
    #filter out lost_sales
    due_dict = Dict(orders_df.id .=> orders_df.due)
    @chain fulfillments_df begin
        @rsubset(
            :delivered <= due_dict[:id], #remove any late fulfillments
            :sent != Symbol("lost_sale"), #remove any lost sales
            view = true
        )
        groupby([:id,:material,:arc])
        @combine(:fulfilled = sum(:amount)) #sum deliveries on each order
        leftjoin(orders_df, _, on = [:id,:material,:arc]) #NOTE: will exclude reallocated orders (have a different arc on the fulfillments df)
        @rtransform(
            :fill = :fulfilled / :amount, #percentage filled
            :service = :fulfilled == :amount #fulfilled?
        )
        coalesce.(0) #replace missings
        groupby([:material,:arc])
        @combine(
            :avg_service = mean(:service),
            :avg_fill = mean(:fill)
        )
    end
end
function calculate_service_measures!(env::SupplyChainEnv; window::Tuple)
    if window == (0,Inf)
        orders_df = env.orders
        fulfillments_df = env.fulfillments
    else
        orders_df = filter(:created => t -> window[1] <= t <= window[2], env.orders, view = true)
        fulfillments_df = filter(:id => in(Set(orders_df.id)), env.fulfillments, view=true)
    end
    env.metrics = calculate_service_measures(orders_df, fulfillments_df)
end