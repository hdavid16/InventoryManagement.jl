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
    check_inputs!(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int)

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs!(network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
                    mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int)
    #initialize flags
    truncate_flag = false
    roundoff_flag1 = false
    roundoff_flag2 = false
    replace_flag = false

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
    
    #run checks on the parameter inputs
    nonsources = [n for n in nodes if !isempty(inneighbors(network, n))]
    env_keys = map_env_keys(nodes, arcs, mrkts, plants, nonsources)
    for obj in keys(env_keys), key in env_keys[obj]
        !in(key, keys(props(network, obj...))) && set_prop!(network, obj..., key, Dict()) #create empty params for object if not specified
        param_dict = get_prop(network, obj..., key)
        for mat in mats
            #set defaults
            if !in(mat, keys(param_dict)) #if material not specified, add it to the dict and set default values
                if key in [:inventory_capacity, :production_capacity] #default is uncapacitated
                    set_prop!(network, obj, key, merge(param_dict, Dict(mat => Inf)))
                elseif key in [:demand_distribution, :lead_time] #zero demand/lead_time for that material
                    set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [0])))
                elseif key == :demand_sequence #zeros demand sequence for that material
                    merge!(param_dict, Dict(mat => zeros(num_periods)))
                elseif key == :demand_frequency #demand at every period
                    merge!(param_dict, Dict(mat => 1))
                elseif key == :supplier_priority #random ordering of supplier priority
                    merge!(param_dict, Dict(mat => inneighbors(network, obj)))
                else #all others default to 0
                    merge!(param_dict, Dict(mat => 0.))
                end
            end
            #check parameter value
            param_dict = get_prop(network, obj..., key)
            param = param_dict[mat]
            if key in [:demand_distribution, :lead_time]
                #convert deterministic value to singleton
                if param isa Number 
                    set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [param])))
                end
                param_dict = get_prop(network, obj..., key)
                param = param_dict[mat]
                #checks
                @assert mean(param) >= 0 "Parameter $key for material $mat at $obj cannot have a negative mean."
                @assert rand(param) isa Number "Parameter $key for material $mat at $obj must be a sampleable distribution or an array."
                param isa Array && @assert length(param) == 1 "Parameter $key for material $mat at $obj cannot be an Array with more than 1 element."
                #convert to deterministic if no variance or truncate if distribution is negative-valued
                if iszero(std(param)) #will only catch Distribution with no std; std(Number) = NaN
                    set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [mean(param)])))
                    replace_flag = true
                elseif minimum(param) < 0 
                    set_prop!(network, obj..., key, merge(param_dict, Dict(mat => truncated(param, 0, Inf))))
                    truncate_flag = true
                end
                if key == :lead_time && (param isa Distribution{Univariate, Continuous} || (param isa Array && mod(param[1],1) > 0))
                    roundoff_flag2 = true
                end    
            elseif key == :demand_sequence
                @assert length(param) == num_periods "The demand sequence for material $mat at node $obj must be a vector with $num_periods entries."
            elseif key == :production_time && mod(param,1) > 0 #check for round-off error in production time
                roundoff_flag1 = true  
            elseif key == :supplier_priority
                for s in param
                    @assert s in inneighbors(network, obj) "Supplier $s is not a supplier to node $obj, but is listed in the supplier priority for that node for material $mat."
                end
            else
                @assert param >= 0 "Parameter $key for material $mat at $obj must be non-negative."
            end
        end
    end

    # for n in nodes, key in all_keys
    #     !in(key, keys(props(network, n))) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
    #     param_dict = get_prop(network, n, key)
    #     for p in mats
    #         if !in(p, keys(param_dict)) #if material not specified, add it to the dict and set default values
    #             if key == :inventory_capacity #default is uncapacitated inventory
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => Inf)))
    #             else #others default to zero
    #                 merge!(param_dict, Dict(p => 0))
    #             end
    #         end
    #         param = param_dict[p]
    #         @assert param >= 0 "Parameter $key for material $p at node $n must be non-negative."
    #     end
    # end
    
    # for n in mrkts, key in market_keys
    #     !in(key, keys(props(network, n))) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
    #     param_dict = get_prop(network, n, key)
    #     for p in mats
    #         if !in(p, keys(param_dict)) #if material not specified, add it to the dict and set default values
    #             if key == :demand_distribution #zero demand for that material
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => [0])))
    #             elseif key == :demand_sequence #zeros demand sequence for that material
    #                 merge!(param_dict, Dict(p => zeros(num_periods)))
    #             elseif key == :demand_frequency #demand at every period (default)
    #                 merge!(param_dict, Dict(p => 1))
    #             else #others default to 0
    #                 merge!(param_dict, Dict(p => 0))
    #             end
    #         end
    #         param = param_dict[p]
    #         if key == :demand_distribution
    #             if param isa Number #convert deterministic value to singleton
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => [param])))
    #             end
    #             @assert mean(param) >= 0 "Parameter $key for material $p at node $n cannot have a negative mean."
    #             @assert rand(param) isa Number "Parameter $key for material $p at node $n must be a sampleable distribution or an array."
    #             param isa Array && @assert length(param) == 1 "Parameter $key for material $p at node $n cannot be an Array with more than 1 element."
    #             if iszero(std(param)) #will only catch Distribution with no std; std(Number) = NaN
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => [mean(param)])))
    #                 replace_flag = true
    #             elseif minimum(param) < 0 
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => truncated(param, 0, Inf))))
    #                 truncate_flag = true
    #             end
    #         elseif key == :demand_sequence
    #             @assert length(param) == num_periods "The demand sequence for material $p at node $n must be a vector with $num_periods entries."
    #         else
    #             @assert param >= 0 "Parameter $key for material $p at node $n must be non-negative."
    #         end
    #     end
    # end
    
    # for n in plants, key in plant_keys
    #     !in(key, keys(props(network, n))) && set_prop!(network, n, key, Dict()) #create empty params for nodes if not specified
    #     param_dict = get_prop(network, n, key)
    #     for p in mats
    #         if !in(p, keys(param_dict)) #if material not specified, add it to the dict and set its default value
    #             if key == :production_capacity #default is uncapacitated production
    #                 set_prop!(network, n, key, merge(param_dict, Dict(p => Inf)))
    #             else #others default to zero
    #                 merge!(param_dict, Dict(p => 0))
    #             end
    #         end
    #         param = param_dict[p]
    #         @assert param >= 0 "Parameter $key for material $p at node $n must be non-negative."
    #         if key == :production_time && mod(param,1) > 0
    #             roundoff_flag1 = true
    #         end
    #     end
    # end
    
    # for a in arcs, key in arc_keys
    #     !in(key, keys(props(network, a...))) && set_prop!(network, a..., key, Dict()) #create empty params for arcs if not specified
    #     param_dict = get_prop(network, a..., key)
    #     for p in mats
    #         if !in(p, keys(param_dict)) #if material not specified, add it to the dict and set its value to 0
    #             if key == :lead_time
    #                 set_prop!(network, a..., key, merge(param_dict, Dict(p => [0])))
    #             else
    #                 merge!(param_dict, Dict(p => 0.))
    #             end
    #         end     
    #         param = param_dict[p]       
    #         if key == :lead_time
    #             param = param_dict[p]
    #             if param isa Number #convert deterministic value to singleton
    #                 set_prop!(network, a..., key, merge(param_dict, Dict(p => [param])))
    #             end
    #             @assert mean(param) >= 0 "Parameter $key for material $p on arc $a cannot have a negative mean."
    #             @assert rand(param) isa Number "Parameter $key for material $p on arc $a must be a sampleable distribution or an array."
    #             param isa Array && @assert length(param) == 1 "Parameter $key for material $p on arc $a cannot be an Array with more than 1 element."
    #             if iszero(std(param)) #convert Distribution with no variance to a deterministic singleton with its mean (will only catch Distribution since std of a singleton array is NaN)
    #                 set_prop!(network, a..., key, merge(param_dict, Dict(p => [mean(param)])))
    #                 replace_flag = true
    #             elseif minimum(param) < 0 
    #                 set_prop!(network, a..., key, merge(param_dict, Dict(p => truncated(param, 0, Inf))))
    #                 truncate_flag = true
    #             end
    #             if param isa Distribution{Univariate, Continuous} || (param isa Array && mod(param[1],1) > 0)
    #                 roundoff_flag2 = true
    #             end    
    #         else
    #             @assert param >= 0 "Parameter $key for material $p on arc $a must be non-negative."
    #         end
    #     end
    # end
    
    # nonsources = [n for n in nodes if !isempty(inneighbors(network, n))]
    # for n in nonsources, key in [:supplier_priority]
    #     !in(key, keys(props(network, n))) && set_prop!(network, n, key, Dict()) #create empty params for market nodes if not specified
    #     param_dict = get_prop(network, n, key)
    #     for p in mats
    #         if !in(p, keys(param_dict)) #if material not specified, add it to the dict and add the node's predecessors
    #             merge!(param_dict, Dict(p => inneighbors(network, n)))
    #         end
    #         param = param_dict[p]
    #         for s in param
    #             @assert s in inneighbors(network, n) "Supplier $s is not a supplier to node $n, but is listed in the supplier priority for that node for material $p."
    #         end
    #     end
    # end

    truncate_flag && @warn "One or more probabilistic distributions allows negative values. The distribution(s) will be truncated to allow only positive values. Note, this will shift the mean of the truncated distribution(s)."
    roundoff_flag1 && @warn "One or more production times are not integer. Round-off error will occur because the simulation uses discrete time."
    roundoff_flag2 && @warn "One or more lead time distributions are not discrete. Round-off error will occur because the simulation uses discrete time."
    replace_flag && @warn "One or more probabilistic distributions passed has zero variance. It has been replaced with its mean value."
end

"""
    map_env_keys(nodes::Vector, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector)

Create a dictionary with the parameter keys for each node/arc in the network
"""
function map_env_keys(nodes::Base.OneTo, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector)
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost]
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :demand_penalty, :demand_sequence]
    plant_keys = [:production_cost, :production_time, :production_capacity]
    arc_keys = [:sales_price, :transportation_cost, :pipeline_holding_cost, :lead_time]
    all_market_keys = vcat(all_keys, market_keys)
    all_plant_keys = vcat(all_keys, plant_keys)

    #list of nodes and arcs
    env_obj = vcat(nodes, arcs)
    #assign keys to each object
    env_keys = Dict(
        obj => 
            obj in nodes ? 
                obj in mrkts ? all_market_keys : 
                obj in plants ? all_plant_keys : 
                all_keys : 
            arc_keys
        for obj in env_obj
    )
    #add supplier priority to the keys for the nonsources
    for node in nodes
        if node in nonsources
            push!(env_keys[node], :supplier_priority)
        end
    end

    return env_keys
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
    for mat0 in mats, mat in setdiff(mats, [mat0])
        i = findfirst(x -> x == mat0, mats)
        j = findfirst(x -> x == mat, mats)
        yen_k = yen_k_shortest_paths(mat_graph, i, j, weights(mat_graph), 100)
        stoich = 0
        for arr in yen_k.paths
            stoich += prod([get_prop(mat_graph, arr[k], arr[k+1], :weight) for k in 1:length(arr)-1])
        end
        mat_conv[(mat0, mat)] = -abs(stoich)
    end

    return mat_graph, mat_conv
end