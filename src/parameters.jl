"""
    identify_nodes(net::MetaDiGraph)

Identify market and producer nodes in the network.
"""
function identify_nodes(net::MetaDiGraph)
    nodes = vertices(net)
    #get end distributors, producers, and distribution centers
    market_keys = [:demand_distribution, :demand_period, :sales_price, :stockout_penalty, :demand_sequence] #keys to identify a market
    plant_keys = [:bill_of_materials, :production_capacity] #keys to identify a plant (producer)
    mrkts = [n for n in nodes if !isempty(intersect(market_keys, keys(net.vprops[n])))]
    plants = [n for n in nodes if !isempty(intersect(plant_keys, keys(net.vprops[n])))]

    return mrkts, plants
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
    check_inputs!(
        network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
        mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int
    )

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs!(
    network::MetaDiGraph, nodes::Base.OneTo, arcs::Vector,
    mrkts::Vector, plants::Vector, mats::Vector, num_periods::Int
)
    #initialize flags
    truncate_flag, roundoff_flag, replace_flag = false, false, false
    
    #run checks on the parameter inputs
    nonsources = [n for n in nodes if !isempty(inneighbors(network, n))] #list of nodes with predecessors
    env_keys = map_env_keys(nodes, arcs, mrkts, plants, nonsources) #dictionary with parameter keys for each node
    for obj in keys(env_keys), key in env_keys[obj]
        if key == :bill_of_materials #check bill of materials (node specific)
            check_bill_of_materials!(network, obj)
        else
            !in(key, keys(props(network, obj...))) && set_prop!(network, obj..., key, Dict()) #create empty params for object if not specified
            param_dict = get_prop(network, obj..., key) #parameter dictionary
            for mat in mats
                #set defaults
                if !in(mat, keys(param_dict)) #if material not specified, add it to the dict and set default values
                    param_dict = set_default!(network, key, obj, mat, num_periods)
                end
                #check parameter value is valid
                param = param_dict[mat]
                if key == :supplier_priority
                    param_dict = check_supplier_priority(network, obj, mat)
                elseif key in [:demand_distribution, :lead_time]
                    param_dict, replace_flag, truncate_flag, roundoff_flag = check_stochastic_variable!(network, key, obj, mat)
                elseif key == :demand_sequence
                    @assert length(param) == num_periods "The demand sequence for material $mat at node $obj must be a vector with $num_periods entries."
                else
                    @assert param isa Real && param >= 0 "Parameter $key for material $mat at $obj must be a non-negative `Real`."
                end
            end
        end
    end

    truncate_flag && @warn """
    One or more probabilistic distributions allows negative values. 
    The distribution(s) will be truncated to allow only positive values. 
    Note, this will shift the mean of the truncated distribution(s).
    """
    roundoff_flag && @warn """
    One or more lead time distributions are not discrete. 
    Round-off error will occur because the simulation uses discrete time.
    """
    replace_flag && @warn """
    One or more probabilistic distributions passed has zero variance. 
    It has been replaced with its mean value.
    """
end

"""
    map_env_keys(nodes::Vector, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector)

Create a dictionary with the parameter keys for each node/arc in the network
"""
function map_env_keys(nodes::Base.OneTo, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector)
    #lists of parameter keys
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost, :service_time]
    market_keys = [:demand_distribution, :demand_period, :sales_price, :stockout_penalty, :demand_sequence]
    plant_keys = [:bill_of_materials, :production_capacity]
    arc_keys = [:sales_price, :stockout_penalty, :transportation_cost, :pipeline_holding_cost, :lead_time]
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
    check_bill_of_materials!(network::MetaDiGraph, n::Int)

Validate bill of material input at a node.
"""
function check_bill_of_materials!(network::MetaDiGraph, n::Int)
    param = get_prop(network, n, :bill_of_materials)
    #create a NamedArray from a BOM Dictionary
    if param isa Dict
        i = unique(first.(keys(param))) #input materials
        j = unique(last.(keys(param))) #output materials
        bom = NamedArray( #initialize named array
            zeros(length(i),length(j)),
            (i,j),
            (:in,:out)
        )
        for key in keys(param) #populate array
            bom[key...] = param[key]
        end
        set_prop!(network, n, :bill_of_materials, bom) #update network metadata
        param = get_prop(network, n, :bill_of_materials)
    end
    #validate
    mats = get_prop(network, :materials)
    @assert (names(param,1) ⊆ mats) && (names(param,2) ⊆ mats) "Bill of material at node $n contains material names that have not been specified in the network metadata."
    @assert param isa NamedArray "The bill of materials at node $n must be a NamedArray."
end

"""
    set_default!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Symbol, num_periods::Int)

Set parameter defaults.
"""
function set_default!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Symbol, num_periods::Int)
    param_dict = get_prop(network, obj..., key)
    if key in [:inventory_capacity, :production_capacity] #default is uncapacitated
        set_prop!(network, obj, key, merge(param_dict, Dict(mat => Inf)))
    elseif key in [:demand_distribution, :lead_time] #zero demand/lead_time for that material
        set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [0])))
    elseif key == :demand_sequence #zeros demand sequence for that material
        merge!(param_dict, Dict(mat => zeros(num_periods)))
    elseif key == :demand_period #demand at every period
        merge!(param_dict, Dict(mat => 1))
    elseif key == :supplier_priority #random ordering of supplier priority
        merge!(param_dict, Dict(mat => inneighbors(network, obj)))
    else #all others default to 0
        merge!(param_dict, Dict(mat => 0.))
    end
    param_dict = get_prop(network, obj..., key)

    return param_dict
end

"""
    check_supplier_priority(network::MetaDiGraph, obj::Int, mat::Symbol)

Validate supplier priority input.
"""
function check_supplier_priority(network::MetaDiGraph, obj::Int, mat::Symbol)
    param_dict = get_prop(network, obj, :supplier_priority)
    param = param_dict[mat]
    #if singleton is provided, put it in a Vector
    if param isa Number 
        set_prop!(network, obj, :supplier_priority, merge(param_dict, Dict(mat => [param])))
    end
    #check suppliers are connected to that node
    param_dict = get_prop(network, obj..., :supplier_priority)
    @assert param_dict[mat] ⊆ inneighbors(network, obj) "One or more suppliers listed in the supplier priority of node $obj is not a predecessor of that node."

    return param_dict
end

"""
    check_stochastic_variable!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Symbol)

Validate inputs for stochastic variables (lead time, demand quantity).
"""
function check_stochastic_variable!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Symbol)
    param_dict = get_prop(network, obj..., key)
    param = param_dict[mat]
    #if singleton is provided, put it in a Vector
    if param isa Number 
        set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [param])))
    end
    param_dict = get_prop(network, obj..., key) #update
    param = param_dict[mat]
    #checks
    @assert param isa Sampleable || (param isa Vector && length(param) == 1 && param[1] >= 0) """
        Parameter $key for material $mat at $obj must be a Sampleable distribution 
        or a sigleton (non-negative) Vector.
    """
    #convert to deterministic if distribution has no variance, or truncate if distribution is negative-valued
    replace_flag, truncate_flag, roundoff_flag = false, false, false
    if iszero(std(param)) #will only catch Distribution with no std; std(Number) = NaN
        set_prop!(network, obj..., key, merge(param_dict, Dict(mat => [mean(param)])))
        replace_flag = true
    elseif minimum(param) < 0 
        set_prop!(network, obj..., key, merge(param_dict, Dict(mat => truncated(param, 0, Inf))))
        truncate_flag = true
    end
    #lead time roundoff will occur if distribution is continuous or deterministic value is not integer
    roundoff_flag = (key == :lead_time) && (param isa Distribution{T, Continuous} where T || (param isa Vector && !iszero(mod.(param,1))))
    #update param dict
    param_dict = get_prop(network, obj..., key) 

    return param_dict, replace_flag, truncate_flag, roundoff_flag
end
