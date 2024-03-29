"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).

[Bug] NOTE: 0 lead time arrivals are not included in `arrivals` (These are updated when the shipment is created)
"""
function calculate_profit!(x::SupplyChainEnv)
    #filter data
    type_set = Set([:delivered, :lost_sale])
    sales_df1 = subset(x.fulfillments, #delivered orders + lost sales
        :period => ByRow(==(x.period)),
        :type => ByRow(in(type_set)),
        view=true
    )
    id_set = Set(sales_df1.id)
    orders_df1 = subset(x.open_orders, #backlogged order quantities
        :id => ByRow(in(id_set)),
        :due => ByRow(<=(0)),
        view=true
    )
    sales_grp1 = groupby(sales_df1, [:arc, :material])
    orders_grp1 = groupby(orders_df1, [:arc, :material])

    #evaluate node profit
    for n in vertices(x.network)
        profit = 0. #initialize node profit
        for mat in x.materials
            #holding cost
            profit += holding_costs(x, n, mat)
            #sales profit at markets (and penalize for unfulfilled demand)
            if n in x.markets
                profit += sum(sales_and_penalties(x, n, :market, mat, sales_grp1, orders_grp1))
            end
            #pay suppliers for received inventory and pay transportation/production cost
            for sup in inneighbors(x.network, n)
                profit += purchases_and_pipeline_costs(x, sup, n, mat, sales_grp1)
            end
            #receive payment for delivered inventory and pay unfulfilled demand penalties
            for req in outneighbors(x.network, n)
                profit += sum(sales_and_penalties(x, n, req, mat, sales_grp1, orders_grp1))
            end
        end

        #discounted profit
        push!(x.profit, [x.period, 1/(1+x.discount)^x.period * profit, n])
    end
end

"""
    holding_costs(x::SupplyChainEnv, n::Int, mat::Material)

Calculate inventory holding costs (negative).
"""
function holding_costs(x::SupplyChainEnv, n::Int, mat::Material)
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
    sales_and_penalties(x::SupplyChainEnv, n::Int, req::Union{Int,Symbol}, mat::Material, sales_grp1::GroupedDataFrame, orders_grp1::GroupedDataFrame)

Calculate sales (positive) and unfulfilled demand penalties (negative) for demand.
"""
function sales_and_penalties(x::SupplyChainEnv, n::Int, req::Union{Int,Symbol}, mat::Material, sales_grp1::GroupedDataFrame, orders_grp1::GroupedDataFrame)
    s, p = 0, 0 #sales and penalties
    arc = req == :market ? (n,) : (n,req)
    sales_price = get_prop(x.network, arc..., :sales_price)[mat]
    dmnd_penalty = get_prop(x.network, arc..., :unfulfilled_penalty)[mat]
    key = (arc = (n,req), material = mat)
    if (sales_price > 0 || dmnd_penalty > 0) && key in keys(sales_grp1)
        sold = sum(filter(:type => !=(:lost_sale), sales_grp1[key], view=true).amount; init=0)
        if key in keys(orders_grp1)
            unfilled = orders_grp1[key].amount[1]
        else
            unfilled = sum(filter(:type => ==(:lost_sale), sales_grp1[key], view=true).amount; init=0)
        end
        s += sales_price * sold
        p -= dmnd_penalty * unfilled
    end

    return s, p
end

"""
    purchases_and_pipeline_costs(x::SupplyChainEnv, sup::Int, n::Int, mat::Material, sales_grp1::GroupedDataFrame)

Calculate payments to suppliers and transportation costs (fixed and variable) (negative).
"""
function purchases_and_pipeline_costs(x::SupplyChainEnv, sup::Int, n::Int, mat::Material, sales_grp1::GroupedDataFrame)
    price = get_prop(x.network, sup, n, :sales_price)[mat]
    trans_cost = get_prop(x.network, sup, n, :transportation_cost)[mat]
    pipe_holding_cost = get_prop(x.network, sup, n, :pipeline_holding_cost)[mat]
    ap = 0 #accounts payable
    key = (arc = (sup,n), material = mat)
    if (price > 0 || trans_cost > 0) && key in keys(sales_grp1) #pay purchase of inventory and transportation cost (assume it is paid to a third party)
        purchased = sum(filter(:type => !=(:lost_sale), sales_grp1[key], view=true).amount; init=0)
        ap -= purchased * price
        ap -= purchased * trans_cost
    end
    if pipe_holding_cost > 0 #pay pipeline holding cost (paid for in-transit inventory in the current period)
        intransit = x.tmp[(sup,n),mat,:pipeline]
        ap -= intransit * pipe_holding_cost
    end

    return ap
end