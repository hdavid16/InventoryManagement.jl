"""
    show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)

Visualize a replenishment order vector as a DataFrame indicating how much of
which material is being requested on each arc.
"""
function show_action(action::Vector{T} where T <: Real, env::SupplyChainEnv)
    mats = env.materials
    arcs = [(e.src, e.dst) for e in edges(env.network)]
    act = reshape(action, (length(mats), length(arcs)))
    return DataFrame(:material => mats, 
                     [Symbol(a) => act[:,i] for (i,a) in enumerate(arcs)]...)
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
    arcs = [(e.src, e.dst) for e in edges(x.network)]
    act = NamedArray(
        reshape(action, (length(x.materials), length(arcs))),
        (x.materials, arcs),
        (:material, :arc)
    )
    #increase period counter
    x.period += 1
    #intialize next period on-hand and pipeline inventories with previous inventories
    initialize_inventories!(x)
    #move active shipments forward one period
    x.shipments.lead .-= 1
    #decrease time until due date for active orders by one period
    x.orders.due .-= 1
    #place requests
    place_orders!(x, act, arcs)
    #update on hand and pipeline inventories due to arrived shipments
    arrivals = update_shipments!(x)
    #discard any excess inventory
    x.options[:capacitated_inventory] && enforce_inventory_limits!(x)
    #markets open and demand occurs
    simulate_markets!(x)
    #update inventories at each node and echelon
    update_inventories!(x)
    #calculate profit at each node
    if x.options[:evaluate_profit]
        calculate_profit!(x, arrivals)
        x.reward = sum(filter(:period => j -> j == x.period, x.profit, view=true).value) #update reward (current profit). NOTE: is this ok for RL?
    end
end

"""
    initialize_inventories!(x::SupplyChainEnv)

Initialize on-hand and pipeline inventories with those from the previous period.
"""
function initialize_inventories!(x::SupplyChainEnv)
    #find previous period's inventories
    prev_inv_on_hand = filter(:period => j -> j == x.period-1, x.inv_on_hand, view=true) #previous inventory level
    prev_onhand_grp = groupby(prev_inv_on_hand, [:node, :material])
    prev_inv_pipeline = filter(:period => j -> j == x.period-1, x.inv_pipeline, view=true) #previous inventory level
    prev_pipeln_grp = groupby(prev_inv_pipeline, [:arc, :material])
    #update dataframes
    for mat in x.materials
        for n in vertices(x.network)
            prev = prev_onhand_grp[(node = n, material = mat)][1,:level] #previous on hand inventory
            push!(x.inv_on_hand, [x.period, n, mat, prev, 0]) #intialize with previous inventory levels
        end
        for a in edges(x.network)
            prev = prev_pipeln_grp[(arc = (a.src,a.dst), material = mat)][1,:level] #previous pipeline inventory
            push!(x.inv_pipeline, [x.period, (a.src,a.dst), mat, prev]) #intialize with previous inventory levels
        end
    end
end

"""
    place_orders!(x::SupplyChainEnv, act::NamedArray, arcs::Vector)

Place inventory replenishment requests throughout the network.
"""
function place_orders!(x::SupplyChainEnv, act::NamedArray, arcs::Vector)
    #exit if no action
    if iszero(act)
        exit_order!(x, arcs)
        return
    end
    #extract info
    mats = x.materials #materials
    #store original production capacities (to account for commited capacity and commited inventory in next section)
    capacities = Dict(n => get_prop(x.network, n, :production_capacity) for n in x.producers)
    #sample lead times
    leads = Dict((a,mat) => rand(get_prop(x.network, a, :lead_time)[mat]) for a in edges(x.network), mat in mats)
    #non source nodes
    nonsources = [n for n in vertices(x.network) if !isempty(inneighbors(x.network, n))]
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #place requests
    for mat in mats #loop by materials
        for req in nonsources #loop by nodes placing requests
            sup_priority = get_prop(x.network, req, :supplier_priority)[mat] #get supplier priority list
            for sup in sup_priority #loop by supplier priority
                a = (sup, req) #arc
                lead = float(leads[Edge(a...), mat]) #sampled lead time
                amount = act[mat, a] #amount requested
                #calculate backlog (will be 0 if assuming lost sales)
                backlog = calculate_backlog(x, sup, req, mat)
                #continue to next iteration if no request made and no backlog
                if iszero(amount) && iszero(backlog)
                    push!(x.demand, [x.period, a, mat, 0, 0, 0, 0, missing])
                    continue 
                end
                #create order and save service time
                x.num_orders += 1 #create new order
                push!(x.all_orders, [x.num_orders, x.period, a, mat, amount, []]) #update order history
                service_time = get_prop(x.network, sup, :service_time)[mat] #get service time
                push!(x.orders, [x.num_orders, a, mat, amount, service_time]) #add order to temp order df
                sort!(x.orders, [:due, :id]) #sort by service time and creation date (serve orders with lowest service time, ranked by order age)
                #try to fulfill order from stock
                fulfill_from_stock!(x, a..., mat, lead, supply_grp, pipeline_grp)
                #try to fulfill order from production
                if sup in x.producers
                    fulfill_from_production!(x, a..., mat, lead, supply_grp, pipeline_grp, capacities)
                end
            end
        end
    end

    #updated production capacities (after commited production)
    for n in x.producers
        set_prop!(x.network, n, :production_capacity, capacities[n])
    end
end

"""
    exit_order!(x, arcs)

Abort order placement.
"""
function exit_order!(x::SupplyChainEnv, arcs::Vector)
    for a in arcs, mat in x.materials
        backlog = calculate_backlog(x, a..., mat)
        push!(x.demand, [x.period, a, mat, 0, 0, 0, backlog, missing])
    end
end

"""
    calculate_backlog(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol)

Calculate quantity of material `mat` backlogged on the arc (`sup`,`req`).
"""
function calculate_backlog(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol)
    due_orders = filter([:arc, :material, :due] => (a,m,t) -> t <= 0 && a == (sup,req) && m == mat, x.orders, view=true)
    backlog = reduce(+, due_orders.quantity; init = 0)
    
    return backlog
end

"""
    fulfill_from_stock!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame
    )

Fulfill request from on-hand inventory.
"""
function fulfill_from_stock!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame
)
    #available supply
    supply = supply_grp[(node = src, material = mat)].level[1]
    
    #loop through orders
    orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.orders, view = true)
    for row in eachrow(orders_df)
        order_amount = row.quantity
        accepted_inv = min(order_amount, supply) #amount fulfilled from inventory
        if accepted_inv > 0
            row.quantity -= accepted_inv #update x.orders (deduct fulfilled part of the order)
            push!(x.all_orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_inv)) #update fulfilled column in order history (date, supplier, amount fulfilled)
            supply_grp[(node = src, material = mat)].level[1] -= accepted_inv #remove inventory from site
            make_shipment!(x, src, dst, mat, accepted_inv, order_amount, lead, pipeline_grp) #ship material
        end
        if row.quantity > 0 && row.due <= 0 && !in(src, x.producers) #if some amount of the order is due and wasn't fulfilled try to reallocate (only if it is not a plant since the plant will try to fulfill from production next)
            update_demand_df!(x, src, dst, mat, accepted_inv, row)
        end
    end
    #remove any fulfilled orders from x.orders
    filter!(:quantity => q -> q > 0, x.orders) 
end

"""
    fulfill_from_production!(
        x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
        lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
    )

Fulfill request by scheduling material production.
"""
function fulfill_from_production!(
    x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, 
    lead::Float64, supply_grp::GroupedDataFrame, pipeline_grp::GroupedDataFrame, capacities::Dict
)
    #extract info
    bom = get_prop(x.network, src, :bill_of_materials)
    rmats = findall(k -> k < 0, bom[:,mat]) #indices of raw materials involved with production of mat
    cmats = findall(k -> k > 0, bom[:,mat]) #indices for co-products
    rmat_names = names(bom,1)[rmats] #names of raw materials
    cmat_names = names(bom,1)[cmats] #names of co-products

    #loop through orders
    orders_df = filter([:arc, :material] => (a,m) -> a == (src,dst) && m == mat, x.orders, view = true)
    for row in eachrow(orders_df)
        cap_and_sup = get_capacity_and_supply(src, mat, bom, rmat_names, cmat_names, capacities, supply_grp)
        order_amount = row.quantity #amount requested in order
        accepted_prod = min(order_amount, cap_and_sup...) #fulfill remaining with available capacity
        if accepted_prod > 0
            row.quantity -= accepted_prod #update x.orders (deduct fulfilled quantity)
            push!(x.all_orders[row.id,:fulfilled], (time=x.period, supplier=src, amount=accepted_prod)) #update fulfilled column in order history (date, supplier, and amount fulfilled)
            capacities[src][mat] -= accepted_prod #update production capacity to account for commited capacity (handled first come first serve)
            for rmat in rmat_names #reactant consumed
                supply_grp[(node = src, material = rmat)].level[1] += accepted_prod * bom[rmat,mat]
            end
            make_shipment!(x, src, dst, mat, accepted_prod, order_amount, lead, pipeline_grp) #schedule shipment
            for cmat in cmat_names #coproducts scheduled for production
                coproduction = accepted_prod*bom[cmat,mat]
                capacities[src][cmat] -= coproduction #update coproduct production capacity
                make_shipment!(x, src, dst, cmat, coproduction, 0, lead, pipeline_grp) #schedule shipment of co-product (original order amount is 0 since this is a coproduct)
            end
        end
        if row.quantity > 0 && row.due <= 0 #if some amount of the order is due and wasn't fulfilled
            update_demand_df!(x, src, dst, mat, accepted_prod, row)
        end
    end
    #remove any fulfilled orders from x.orders
    filter!(:quantity => q -> q > 0, x.orders) 
end

"""
    get_capacity_and_supply(
        bom::NamedArray, rmat_names::Vector, cmat_names::Vector, 
        capacities::Dict, supply_grp::GroupedDataFrame
    )

Get available capacity and material supply at producer.
"""
function get_capacity_and_supply(
    n::Int, mat::Symbol, bom::NamedArray, 
    rmat_names::Vector, cmat_names::Vector, 
    capacities::Dict, supply_grp::GroupedDataFrame
)
    #commit production at plant
    capacity = [] #get production capacity
    mat_supply = [] #store max capacity based on raw material consumption for each raw material
    isempty(rmat_names) && return [0],[0] #if material is not produced at the node, then return zero capacity and zero supply
    push!(capacity, capacities[n][mat]) #production capacity for that material
    for rmat in rmat_names #check raw material supply
        sup_pp = supply_grp[(node = n, material = rmat)].level[1] #supply of material involved in BOM
        push!(mat_supply, - sup_pp / bom[rmat,mat]) #only account for raw materials that are in the BOM
    end 
    for cmat in cmat_names #add capacity constraint for any co-products (scaled by stoichiometry)
        push!(capacity, capacities[n][rmat] / bom[cmat,mat])
    end

    return vcat(capacity, mat_supply)
end

"""
    make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, original_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)

Schedule shipment of material.
"""
function make_shipment!(x::SupplyChainEnv, src::Int, dst::Int, mat::Symbol, accepted_amount::Float64, original_amount::Float64, lead::Float64, pipeline_grp::GroupedDataFrame)
    push!(x.shipments, [(src,dst), mat, accepted_amount, lead])
    push!(x.demand, [x.period, (src,dst), mat, original_amount, accepted_amount, lead, 0, missing])
    pipeline_grp[(arc = (src, dst), material = mat)].level[1] += accepted_amount #add inventory to pipeline
end

"""
    update_demand_df!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, accepted::Float64, order_row::DataFrameRow)

Reallocate order to next supplier in the priority list.
"""
function update_demand_df!(x::SupplyChainEnv, sup::Int, req::Int, mat::Symbol, accepted::Float64, order_row::DataFrameRow)
    sup_priority = get_prop(x.network, req, :supplier_priority)[mat]
    next_sup = findfirst(k -> k == sup, sup_priority) + 1 #get next in line
    new_alloc = missing
    if x.options[:reallocate] && next_sup <= length(sup_priority)
        new_sup = sup_priority[next_sup] #next supplier in line
        new_alloc = (new_sup, req) #store new arc where reallocated
        order_row.arc = new_alloc #update x.orders to reallocate material
    end
    if accepted > 0 #order was partially fulfilled
        x.demand[end,:unfulfilled] = order_row.quantity
        x.demand[end,:reallocated] = new_alloc
    else #order was not fulfilled
        original_amount = order_row.quantity
        push!(x.demand, [x.period, (sup,req), mat, original_amount, 0, 0, original_amount, new_alloc])
    end
end

"""
    update_shipments!(x::SupplyChainEnv)

Update inventories throughout the network for arrived shipments.
"""
function update_shipments!(x::SupplyChainEnv)
    #filter data
    #get on hand inventory
    supply_df = filter(:period => k -> k == x.period, x.inv_on_hand, view=true) #on_hand inventory supply
    supply_grp = groupby(supply_df, [:node, :material]) #group on hand inventory supply
    #get pipeline inventory
    pipeline_df = filter(:period => k -> k == x.period, x.inv_pipeline, view=true)
    pipeline_grp = groupby(pipeline_df, [:arc, :material])

    #find active shipments with 0 lead time
    arrivals = filter(:lead => i -> i <= 0, x.shipments) 
    for i in 1:nrow(arrivals)
        a, mat, amount = arrivals[i, [:arc, :material, :amount]]
        supply_grp[(node = a[end], material = mat)].level[1] += amount
        pipeline_grp[(arc = a, material = mat)].level[1] -= amount
    end
    filter!(:lead => i -> i > 0, x.shipments) #remove shipments that arrived (and any zero values)
    
    return arrivals
end

"""
    enforce_inventory_limits!(x::SupplyChainEnv)

Discard any excess inventory (exceeding the inventory capacity at each node).
"""
function enforce_inventory_limits!(x::SupplyChainEnv)
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventories
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material for easier lookup
    for n in vertices(x.network)
        node_max_inv = get_prop(x.network, n, :inventory_capacity)
        for mat in x.materials
            max_inv = node_max_inv[mat]
            (iszero(max_inv) || isinf(max_inv)) && continue #skip iteration if node doesn't store that material or has infinite capacity
            onhand = onhand_grp[(node = n, material = mat)].level[1] #onhand inventory
            if onhand > max_inv
                onhand_grp[(node = n, material = mat)].level[1] = max_inv
                onhand_grp[(node = n, material = mat)].discarded[1] += onhand - max_inv
            end
        end
    end
end


"""
    update_inventories!(x::SupplyChainEnv)

Update inventory position and inventory level for all materials and nodes.
"""
function update_inventories!(x::SupplyChainEnv)
    #filter data
    on_hand_df = filter(:period => j -> j == x.period, x.inv_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
    pipeline_grp = groupby(pipeline_df, :material) #group by material
    orders_df = filter([:period, :reallocated] => (i1, i2) -> i1 == x.period && ismissing(i2), x.demand, view=true) #replenishment orders from current period
    orders_grp = groupby(orders_df, :material) #group by material

    #initialize echelon inventory positions
    initialize_echelons!(x)
    ech_df = filter(:period => j -> j == x.period, x.ech_position, view=true) #echelon position at each node
    ech_grp = groupby(ech_df, [:node, :material]) #group by material
    
    #loop through nodes and update inventory levels, positions, and echelons
    for n in vertices(x.network), mat in x.materials
        ilevel, ipos = inventory_components(x, n, mat, pipeline_grp, onhand_grp, orders_grp)
        push!(x.inv_level, [x.period, n, mat, ilevel]) #update inventory level
        push!(x.inv_position, [x.period, n, mat, ipos]) #update inventory position
        update_echelons!(x, n, mat, ipos, ech_grp) #update echelon positions
    end

end

"""
    initialize_echelons!(x::SupplyChainEnv)

Initialize echelon positions at the current period to 0.
"""
function initialize_echelons!(x::SupplyChainEnv)
    for n in vertices(x.network), mat in x.materials
        if get_prop(x.network, n, :inventory_capacity)[mat] > 0
            push!(x.ech_position, [x.period, n, mat, 0])
        end
    end
end

"""
    inventory_components(
        x::SupplyChainEnv, n::Int, mat::Symbol, 
        pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, 
        orders_grp::GroupedDataFrame
    )

Extract components to determine inventory level and position.
"""
function inventory_components(
    x::SupplyChainEnv, n::Int, mat::Symbol, 
    pipeline_grp::GroupedDataFrame, onhand_grp::GroupedDataFrame, 
    orders_grp::GroupedDataFrame
)
    pipeline = sum(filter(:arc => j -> j[end] == n, pipeline_grp[(material = mat,)], view=true).level) #in-transit inventory
    onhand = onhand_grp[(node = n, material = mat)].level[1] #on_hand inventory
    backorder = 0 #initialize replenishment orders placed to suppliers that are backlogged
    backlog = 0 #initialize backlog for orders placed by successors
    if x.options[:backlog]
        backlog = sum(filter(:arc => i -> i[1] == n, orders_grp[(material = mat,)], view=true).unfulfilled)
        backorder = sum(filter(:arc => i -> i[2] == n && i[2] != i[1], orders_grp[(material = mat,)], view=true).unfulfilled)
    end
    ilevel = onhand - backlog #inventory level
    iorder = pipeline + backorder #inventory on order
    ipos = ilevel + iorder #inventory position

    return ilevel, ipos
end

"""
    update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ech_grp::GroupedDataFrame)

Update echelon positions for current time period.
"""
function update_echelons!(x::SupplyChainEnv, n::Int, mat::Symbol, ipos::Float64, ech_grp::GroupedDataFrame)
    for ech in findall(i -> n in i, x.echelons) #identify which echelons have been affected and add to these
        if get_prop(x.network, ech, :inventory_capacity)[mat] > 0 #only add to echelon if that node holds that material
            ech_grp[(node = ech, material = mat)].level[1] += ipos
        end
    end
end

"""
    simulate_markets!(x::SupplyChainEnv)

Open markets, apply material demands, and update inventory positions.
"""
function simulate_markets!(x::SupplyChainEnv)
    on_hand_df = filter([:period,:node] => (t,n) -> t == x.period && n in x.markets, x.inv_on_hand, view=true) #on_hand inventory at node
    onhand_grp = groupby(on_hand_df, [:node, :material]) #group by node and material
    last_dmnd = filter([:period,:arc] => (t,a) -> t == x.period-1 && a[1] == a[2], x.demand, view=true)
    last_d_group = groupby(last_dmnd, [:arc, :material])
    for n in x.markets
        dmnd_seq_dict = get_prop(x.network, n, :demand_sequence)
        dmnd_freq_dict = get_prop(x.network, n, :demand_frequency)
        dmnd_dist_dict = get_prop(x.network, n, :demand_distribution)
        for mat in x.materials
            dmnd_seq = dmnd_seq_dict[mat]
            dprob = Bernoulli(1/dmnd_freq_dict[mat]) #demand probability (probability of ordering is 1/demand period; aka, once every x days)
            dmnd = dmnd_dist_dict[mat] #demand distribution
            q = iszero(dmnd_seq) ? rand(dprob) * rand(dmnd) : dmnd_seq[x.period] #quantity requested (sampled or specified by user)
            if x.options[:backlog] && x.period > 1 #add previous backlog to the quantity requested at the market
                q += last_d_group[(arc = (n,n), material = mat)].unfulfilled[1]
            end
            inv = onhand_grp[(node = n, material = mat)].level[1] #on hand inventory
            sold = min(q, inv) #sales
            unfilled = max(q - inv, 0.) #unfilfilled
            push!(x.demand, [x.period, (n,n), mat, q, sold, 0, unfilled, missing]) #update df
            onhand_grp[(node = n, material = mat)].level[1] -= sold #update end of period inventory (subtract sales)
        end
    end
end

"""
    calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)

Calculate profit at each node in the network
    (sales - production cost - purchase costs - transportation costs - holding costs).
"""
function calculate_profit!(x::SupplyChainEnv, arrivals::DataFrame)
    #filter data
    on_hand_df = filter([:period, :level] => (j1, j2) -> j1 == x.period && !isinf(j2), x.inv_on_hand, view=true) #on_hand inventory 
    onhand_grp = groupby(on_hand_df, [:node, :material])
    sales_df = filter([:period, :arc] => (t,a) -> t == x.period && a[1] == a[2], x.demand, view=true) #replenishment orders
    sales_grp = groupby(sales_df, [:arc, :material])
    pipeline_df = filter(:period => j -> j == x.period, x.inv_pipeline, view=true) #pipeline inventories
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