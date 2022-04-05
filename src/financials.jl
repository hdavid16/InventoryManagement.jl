"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    #filter data
    on_hand_df = filter([:period, :level] => (j1, j2) -> j1 == x.period && !isinf(j2), x.inventory_on_hand, view=true) #on_hand inventory 
    onhand_grp = groupby(on_hand_df, [:node, :material])
    sales_df = filter([:period, :arc] => (t,a) -> t == x.period && a[1] == a[2], x.demand, view=true) #replenishment orders
    sales_grp = groupby(sales_df, [:arc, :material])
    pipeline_df = filter(:period => j -> j == x.period, x.inventory_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #evaluate node profit
    for n in vertices(x.network)
        profit = 0. #initialize node profit
        for mat in x.materials
            #holding cost
            profit += holding_costs(x, n, mat, onhand_grp)
            #sales profit at markets (and penalize for unfulfilled demand)
            if n in x.markets
                profit += sum(sales_and_penalties(x, n, mat, sales_grp))
            end
            #pay suppliers for received inventory and pay transportation/production cost
            for pred in inneighbors(x.network, n)
                profit += accounts_payable(x, pred, n, mat, arrivals, pipeline_grp)
            end
            #receive payment for delivered inventory
            for succ in outneighbors(x.network, n)
                profit += accounts_receivable(x, n, succ, mat, arrivals)
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end

"""
    holding_costs(x::SupplyChainEnv, n::Int, mat::Symbol, onhand_grp::GroupedDataFrame)

Calculate holding costs (negative).
"""
function holding_costs(x::SupplyChainEnv, n::Int, mat::Symbol, onhand_grp::GroupedDataFrame)
    #holding cost
    hold_cost = get_prop(x.network, n, :holding_cost)[mat]
    c = 0
    if hold_cost > 0
        onhand = onhand_grp[(node = n, material = mat)].level
        if !isempty(onhand)
            c -= hold_cost * onhand[1]
        end
    end

    return c
end

"""
    sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Symbol, sales_grp::GroupedDataFrame)

Calculate sales (positive) and unfulfilled demand penalties (negative)
"""
function sales_and_penalties(x::SupplyChainEnv, n::Int, mat::Symbol, sales_grp::GroupedDataFrame)
    sales_price = get_prop(x.network, n, :sales_price)[mat]
    dmnd_penalty = get_prop(x.network, n, :demand_penalty)[mat]
    s, p = 0, 0 #sales and penalties
    if sales_price > 0 || dmnd_penalty > 0
        sold, unfilled = sales_grp[(arc = (n,n), material = mat)][1, [:fulfilled, :unfulfilled]]
        s += sales_price * sold
        p -= dmnd_penalty * unfilled
    end

    return s, p
end

"""
    accounts_payable(x::SupplyChainEnv, sup::Int, n::Int, mat::Symbol, arrivals::DataFrame, pipeline_grp::GroupedDataFrame)

Calculate payments to suppliers and transportation costs (negative).
"""
function accounts_payable(x::SupplyChainEnv, sup::Int, n::Int, mat::Symbol, arrivals::DataFrame, pipeline_grp::GroupedDataFrame)
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
        intransit = pipeline_grp[(arc = (sup, n), material = mat)].level[1]
        ap -= intransit * pipe_holding_cost
    end

    return ap
end

"""
    accounts_receivable(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, arrivals::DataFrame)

Calculate payment for sales to downstream nodes (positive).
"""
function accounts_receivable(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, arrivals::DataFrame)
    price = get_prop(x.network, sup, req, :sales_price)[mat]
    ar = 0
    if price > 0 #receive payment for delivered inventory
        sold = filter([:arc, :material] => (j1, j2) -> j1 == (sup,req) && j2 == mat, arrivals, view=true).amount
        if !isempty(sold)
            ar += sold[1] * price
        end
    end

    return ar
end