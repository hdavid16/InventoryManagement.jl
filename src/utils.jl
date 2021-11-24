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
    check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int)

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int)
    truncate_flag = false
    roundoff_flag1 = false
    roundoff_flag2 = false

    #bill of materials
    if :bill_of_materials in keys(network.gprops)
        bom = get_prop(network, :bill_of_materials)
    elseif isempty(plants)
        bom = zeros(Int, length(mats), length(mats))
        set_prop!(network, :bill_of_materials, bom)
    else
        @assert :bill_of_materials in keys(network.gprops) "Bill of materials is missing."
    end
    if bom != [0] #if [0] then single product, no bom
        @assert typeof(bom) <: Matrix{T} where T <: Real "Bill of materials must be a matrix of real numbers."
        @assert size(bom)[1] == size(bom)[2] "Bill of materials must be a square matrix."
    end
    @assert size(bom)[1] == length(mats) "The number of rows and columns in the bill of materials must be equal to the number of materials."
    
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost]
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :demand_penalty, :demand_sequence]
    plant_keys = [:production_cost, :production_time, :production_capacity]
    arc_keys = [:sales_price, :transportation_cost, :pipeline_holding_cost, :lead_time]
    
    for n in nodes, key in all_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set default values
                if key == :inventory_capacity #default is uncapacitated inventory
                    tmp = Dict(p => Inf)
                    network.vprops[n][key] = merge(network.vprops[n][key], tmp)
                else #others default to zero
                    network.vprops[n][key][p] = 0
                end
            end
            @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p at node $n must be non-negative."
        end
    end
    
    for n in mrkts, key in market_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set default values
                if key == :demand_distribution #zero demand for that material
                    tmp = Dict(p => [0])
                    network.vprops[n][key] = merge(network.vprops[n][key], tmp)
                elseif key == :demand_sequence #zeros demand sequence for that material
                    network.vprops[n][key][p] = zeros(num_periods)
                elseif key == :demand_frequency #demand at every period (default)
                    network.vprops[n][key][p] = 1
                else #others default to 0
                    network.vprops[n][key][p] = 0
                end
            end
            if key == :demand_distribution
                dmnd_dst = network.vprops[n][key][p]
                @assert rand(dmnd_dst) isa Number "Parameter $key for material $p at node $n must be a sampleable distribution or an array."
                dmnd_dst isa Array && @assert length(dmnd_dst) == 1 && dmnd_dst[1] >= 0 "Parameter $key for material $p at node $n cannot be negative."
                if minimum(dmnd_dst) < 0 
                    tmp = Dict(p => truncated(dmnd_dst, 0, Inf))
                    network.vprops[n][key] = merge(network.vprops[n][key], tmp)
                    truncate_flag = true
                end
            elseif key == :demand_sequence
                @assert length(network.vprops[n][key][p]) == num_periods "The demand sequence for material $p at node $n must be a vector with $num_periods entries."
            else
                @assert network.vprops[n][key][p] >= 0 "Parameter $key for material $p at node $n must be non-negative."
            end
        end
    end
    
    for n in plants, key in plant_keys
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
        for p in mats
            if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and set its default value
                if key == :production_capacity #default is uncapacitated production
                    network.vprops[n][key][p] = Inf
                else #others default to zero
                    network.vprops[n][key][p] = 0
                end
            end
            param = network.vprops[n][key][p]
            @assert param >= 0 "Parameter $key for material $p at node $n must be non-negative."
            if key == :production_time && mod(param,1) > 0
                roundoff_flag1 = true
            end
        end
    end
    
    for a in arcs, key in arc_keys
        !in(key, keys(network.eprops[Edge(a...)])) && set_prop!(network, Edge(a...), key, Dict()) #create empty params for arcs if not specified
        for p in mats
            if !in(p, keys(network.eprops[Edge(a...)][key])) #if material not specified, add it to the dict and set its value to 0
                if key == :lead_time
                    tmp = Dict(p => [0])
                    network.eprops[Edge(a...)][key] = merge(network.eprops[Edge(a...)][key], tmp)
                else
                    network.eprops[Edge(a...)][key][p] = 0.
                end
            end            
            if key == :lead_time
                lt = network.eprops[Edge(a...)][key][p]
                @assert rand(lt) isa Number "Parameter $key for material $p on arc $a must be a sampleable distribution or an array."
                lt isa Array && @assert length(lt) == 1 && lt[1] >= 0 "Parameter $key for material $p on arc $a cannot be negative and must be a singleton or an univariate distribution."
                if minimum(lt) < 0 
                    tmp = Dict(p => truncated(lt, 0, Inf))
                    network.eprops[Edge(a...)][key] = merge(network.eprops[Edge(a...)][key], tmp)
                    truncate_flag = true
                end
                if lt isa Distribution{Univariate, Continuous} || (lt isa Array && mod(lt[1],1) > 0)
                    roundoff_flag2 = true
                end    
            else
                @assert network.eprops[Edge(a...)][key][p] >= 0 "Parameter $key for material $p on arc $a must be non-negative."
            end
        end
    end
    
    nonsources = [n for n in nodes if !isempty(inneighbors(network, n))]
    for n in nonsources, p in mats
        key = :supplier_priority
        !in(key, keys(network.vprops[n])) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
        if !in(p, keys(network.vprops[n][key])) #if material not specified, add it to the dict and add the node's predecessors
            network.vprops[n][key][p] = inneighbors(network, n)
        end
        for s in network.vprops[n][key][p]
            @assert s in inneighbors(network, n) "Supplier $s is not a supplier to node $n, but is listed in the supplier priority for that node for material $p."
        end
    end

    truncate_flag && @warn "One or more probabilistic distributions allows negative values. The distribution(s) will be truncated to allow only positive values."
    roundoff_flag1 && @warn "One or more production times are not integer. Round-off error will occur because the simulation uses discrete time."
    roundoff_flag2 && @warn "One or more lead time distributions are not discrete. Round-off error will occur because the simulation uses discrete time."
end

"""
    identify_echelons(network::MetaDiGraph, n::Int)

Identify all `network` successors to node `n`.
"""
function identify_echelons(network::MetaDiGraph, n::Int)
    #find all successors to node n
    sink_nodes = [i for i in vertices(network) if isempty(outneighbors(network, i))]
    echelon_nodes = [] #initialize list of echelon nodes
    for sink in sink_nodes #find all nodes between n and each sink node
        echelon = yen_k_shortest_paths(network, n, sink, weights(network), typemax(Int)).paths
        union!(echelon_nodes, echelon...) #remove duplicates
    end

    return echelon_nodes
end

"""
    service_measures(env::SupplyChainEnv)

Calculate mean service level and fill rate for each node and each material in the simulation.

NOTE: Not suported when there is more than 1 supplier and reallocation is allowed.
      CSL is currently on a per-period basis, not a cycle basis
"""
function service_measures(env::SupplyChainEnv; review_period::Union{Int, StepRange, Vector, Dict} = 1)
    #merge demand and replenishment tables
    demand_filt = filter(i -> i.demand > 0, env.demand) #filter out times with no demand at markets
    replenish_filt = filter(i -> i.requested > 0, env.replenishments)
    select!(replenish_filt, #convert replenishment table into same format as demand table for merging
        :period,
        :arc => ByRow(first) => :node,
        :material,
        :requested => :demand,
        :accepted => :sold,
        :unfulfilled
    )
    orders = vcat(demand_filt, replenish_filt) #merge two tables

    #convert orders of products at producer nodes into orders of raw materials at those nodes
    bom = env.bill_of_materials
    for row in eachrow(filter(i -> i.node in env.producers, orders))
        col = findfirst(i -> i == row.material, env.materials)
        for id in findall(i -> i < 0, bom[:,col])
            push!(orders, (row.period, row.node, env.materials[id], -row.demand*bom[id,col], -row.sold*bom[id,col], -row.unfulfilled*bom[id,col]))
        end
    end

    #convert review period into proper format (dictionary by node)
    if review_period isa Int
        review_period = 1:review_period:env.num_periods
    end
    if !isa(review_period, Dict)
        review_period = Dict((n,m) => review_period for n in vertices(env.network), m in env.materials)
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
            [:period,:node,:material] => ByRow((t,n,m) -> findfirst(x -> t in x, cycle[n,m])) => :cycle,
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

"""
    material_conversion(net::MetaDiGraph)

Generate a material graph where the weights are the stoichiometry. 
Also generate a dictionary mapping the amount of each material consumed when one unit of product is procured.

NOTE: assumes that byproducts don't have their own primary pathways
      assumes that byproducts are only made in one of the pathways
"""
function material_conversion(net::MetaDiGraph)
    #get materials
    mats = get_prop(net, :materials)
    num_mats = length(mats) #number of materials
    #expanded bom
    bom = get_prop(net, :bill_of_materials)
    bom_exp = copy(bom) 
    bom_exp[bom .> 0] .= 0 #replace positive numbers (coproduction) with zero
    coprod = findall(i -> i > 0, bom) #find indices where there is coproduction
    for idx in coprod #loop through indices and create columns for coproduction
        coprod_col, copy_col = idx.I
        bom_exp[:, coprod_col] = bom_exp[:, copy_col]
    end
    #generate material graph (and connect nodes)
    mat_graph = MetaDiGraph(num_mats)
    for i in 1:num_mats, j in 1:num_mats
        if bom_exp[i,j] < 0 
            add_edge!(mat_graph, i, j, :weight, bom_exp[i,j])
        end
    end
    #create conversion dictionary
    mat_conv = Dict()
    for m in mats, p in mats
        if m != p
            i = findfirst(x -> x == m, mats)
            j = findfirst(x -> x == p, mats)
            yen_k = yen_k_shortest_paths(mat_graph, i, j, weights(mat_graph), 100)
            if !isempty(yen_k.dists)
                mat_paths = yen_k.paths
                stoich = 0
                for arr in mat_paths
                    stoich += prod([get_prop(mat_graph, arr[k], arr[k+1], :weight) for k in 1:length(arr)-1])
                end
                mat_conv[(m, p)] = -abs(stoich)
            end
        end
    end

    return mat_graph, mat_conv
end