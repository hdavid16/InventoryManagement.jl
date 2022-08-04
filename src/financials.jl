"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).

[Bug] NOTE: 0 lead time arrivals are not included in `arrivals` (These are updated when the shipment is created)
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    #filter data
    sales_df = filter([:period, :arc] => (t,a) -> t == x.period, x.demand, view=true) #replenishment orders
    sales_grp = groupby(sales_df, [:arc, :material])

    #evaluate node profit
    for n in vertices(x.network)
        profit = 0. #initialize node profit
        for mat in x.materials
            #holding cost
            profit += holding_costs(x, n, mat)
            #sales profit at markets (and penalize for unfulfilled demand)
            if n in x.markets
                profit += sum(sales_and_penalties(x, n, mat, sales_grp))
            end
            #pay suppliers for received inventory and pay transportation/production cost
            for sup in inneighbors(x.network, n)
                profit += purchases_and_pipeline_costs(x, sup, n, mat, arrivals)
            end
            #receive payment for delivered inventory and pay unfulfilled demand penalties
            for req in outneighbors(x.network, n)
                profit += sum(sales_and_penalties(x, n, req, mat, sales_grp, arrivals))
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end

"""
    holding_costs(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})

Calculate inventory holding costs (negative).
"""
function holding_costs(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})
    #holding cost
    hold_cost = get_prop(x.network, n, :holding_cost)[mat]
    c = 0
    if hold_cost > 0
        onhand = x.tmp[n,mat,:on_hand]
        if !isinf(onhand)
            c -= hold_cost * onhand[1]
        end
    end

    return c
end

"""
    sales_and_penalties(x::SupplyChainEnv, n::Int, req::Int, mat::Union{Symbol,String}, sales_grp::GroupedDataFrame, arrivals::DataFrame)

Calculate sales (positive) and unfulfilled demand penalties (negative) for internal demand.
"""
function sales_and_penalties(x::SupplyChainEnv, n::Int, req::Int, mat::Union{Symbol,String}, sales_grp::GroupedDataFrame, arrivals::DataFrame)
    s, p = 0, 0 #sales and penalties
    sales_price = get_prop(x.network, n, req, :sales_price)[mat]
    dmnd_penalty = get_prop(x.network, n, req, :unfulfilled_penalty)[mat]
    key = (arc = (n,req), material = mat)
    if (sales_price > 0 || dmnd_penalty > 0) && key in keys(sales_grp)
        delivered = filter([:arc, :material] => (a, m) -> a == (n,req) && m == mat, arrivals, view=true).amount
        sold = !isempty(delivered) ? delivered[1] : 0
        unfilled = sales_grp[key].unfulfilled[1]
        s += sales_price * sold
        p -= dmnd_penalty * unfilled
    end

    return s, p
end

"""
    sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, sales_grp::GroupedDataFrame)

Calculate sales (positive) and unfulfilled demand penalties (negative) for external demand.
"""
function sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String}, sales_grp::GroupedDataFrame)
    s, p = 0, 0 #sales and penalties
    sales_price = get_prop(x.network, n, :sales_price)[mat]
    dmnd_penalty = get_prop(x.network, n, :unfulfilled_penalty)[mat]
    key = (arc = (n,:market), material = mat)
    if (sales_price > 0 || dmnd_penalty > 0) && key in keys(sales_grp)
        sold, unfilled = sales_grp[key][1, [:fulfilled, :unfulfilled]]
        s += sales_price * sold
        p -= dmnd_penalty * unfilled
    end

    return s, p
end

"""
    purchases_and_pipeline_costs(x::SupplyChainEnv, sup::Int, n::Int, mat::Union{Symbol,String}, arrivals::DataFrame)

Calculate payments to suppliers and transportation costs (fixed and variable) (negative).
"""
function purchases_and_pipeline_costs(x::SupplyChainEnv, sup::Int, n::Int, mat::Union{Symbol,String}, arrivals::DataFrame)
    price = get_prop(x.network, sup, n, :sales_price)[mat]
    trans_cost = get_prop(x.network, sup, n, :transportation_cost)[mat]
    pipe_holding_cost = get_prop(x.network, sup, n, :pipeline_holding_cost)[mat]
    ap = 0 #accounts payable
    if price > 0 || trans_cost > 0 #pay purchase of inventory and transportation cost (assume it is paid to a third party)
        purchased = filter([:arc, :material] => (j1, j2) -> j1 == (sup, n) && j2 == mat, arrivals, view=true).amount
        if !isempty(purchased)
            ap -= purchased[1] * price
            ap -= purchased[1] * trans_cost
        end
    end
    if pipe_holding_cost > 0 #pay pipeline holding cost (paid for in-transit inventory in the current period)
        intransit = x.tmp[(sup,n),mat,:pipeline]
        ap -= intransit * pipe_holding_cost
    end

    return ap
end