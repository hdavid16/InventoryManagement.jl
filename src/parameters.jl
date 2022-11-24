"""
    identify_nodes(net::MetaDiGraph)

Identify market and producer nodes in the network.
"""
function identify_nodes(net::MetaDiGraph)
    nodes = vertices(net)
    #get end distributors, producers, and distribution centers
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :unfulfilled_penalty] #keys to identify a market
    plant_keys = [:bill_of_materials, :production_capacity, :make_to_order] #keys to identify a plant (producer)
    mrkts = [n for n in nodes if !isempty(intersect(market_keys, keys(props(net,n))))]
    plants = [n for n in nodes if !isempty(intersect(plant_keys, keys(props(net,n))))]

    return mrkts, plants
end

"""
    identify_echelons(network::MetaDiGraph, n::Int)

Identify all `network` successors to node `n`.
"""
function identify_echelons(network::MetaDiGraph, n::Int)
    net = remove_self_loops(network)
    #find all successors to node n
    sink_nodes = [i for i in vertices(net) if isempty(outneighbors(net, i))]
    echelon_nodes = [] #initialize list of echelon nodes
    for sink in sink_nodes #find all nodes between n and each sink node
        echelon = yen_k_shortest_paths(net, n, sink, weights(net), typemax(Int)).paths
        union!(echelon_nodes, echelon...) #remove duplicates
    end

    return echelon_nodes
end

"""
    check_inputs!(network::MetaDiGraph, mrkts::Vector, plants::Vector)

Check inputs when creating a `SupplyChainEnv`.
"""
function check_inputs!(network::MetaDiGraph, mrkts::Vector, plants::Vector)

    nodes = vertices(network) #network nodes
    arcs = [(e.src,e.dst) for e in edges(network)] #network edges as Tuples (not Edges)
    mats = get_prop(network, :materials) #materials
    sort!(mats)

    #initialize flags
    truncate_flag, roundoff_flag, replace_flag = false, false, false
    
    #run checks on the parameter inputs
    nonsources = filter(n -> !isempty(inneighbors(network, n)), nodes) #list of nodes with predecessors
    nonsinks = filter(n -> !isempty(outneighbors(network, n)), nodes) #list of nodes with successors
    env_keys = map_env_keys(nodes, arcs, mrkts, plants, nonsources, nonsinks) #dictionary with parameter keys for each node
    for obj in keys(env_keys), key in env_keys[obj]
        if key == :bill_of_materials #check bill of materials (node specific)
            check_bill_of_materials!(network, obj)
        elseif key == :make_to_order
            check_make_to_order!(network, obj)
        else
            !(key in keys(props(network, obj...))) && set_prop!(network, obj..., key, Dict()) #create empty params for object if not specified
            param_dict = get_prop(network, obj..., key) #parameter dictionary
            for mat in mats
                #set defaults
                if !(mat in keys(param_dict)) #if material not specified, add it to the dict and set default values
                    param_dict = set_default!(network, key, obj, mat)
                end
                #check parameter value is valid
                param = param_dict[mat]
                if key == :supplier_priority
                    param_dict = check_supplier_priority(network, obj, mat)
                elseif key == :customer_priority
                    param_dict = check_customer_priority(network, obj, mat)
                elseif key in [:demand_distribution, :lead_time, :service_lead_time]
                    param_dict, replace_flag, truncate_flag, roundoff_flag = check_stochastic_variable!(network, key, obj, mat)
                else
                    @assert param isa Real && param >= 0 "Parameter $key for material $mat at $obj must be a non-negative `Real`."
                end
            end
        end
    end

    store_node_materials!(network, plants) #store materials allowed in each node
    
    print_warnings(truncate_flag, roundoff_flag, replace_flag)
end

"""
    map_env_keys(nodes::Vector, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector, nonsinks::Vector)

Create a dictionary with the parameter keys for each node/arc in the network
"""
function map_env_keys(nodes::Base.OneTo, arcs::Vector, mrkts::Vector, plants::Vector, nonsources::Vector, nonsinks::Vector)
    #lists of parameter keys
    all_keys = [:initial_inventory, :inventory_capacity, :holding_cost]
    market_keys = [:demand_distribution, :demand_frequency, :sales_price, :unfulfilled_penalty, :service_lead_time, :market_partial_fulfillment, :market_early_fulfillment]
    plant_keys = [:bill_of_materials, :production_capacity, :make_to_order]
    arc_keys = [:sales_price, :unfulfilled_penalty, :transportation_cost, :pipeline_holding_cost, :lead_time, :service_lead_time]
    all_market_keys = vcat(all_keys, market_keys)
    all_plant_keys = vcat(all_keys, plant_keys)
    all_market_plant_keys = vcat(all_keys, market_keys, plant_keys)
    #list of nodes and arcs
    mrkt_plants = mrkts ∩ plants
    env_obj = vcat(nodes, arcs)
    #assign keys to each object
    env_keys = Dict(
        obj => 
            obj in nodes ? 
                obj in mrkt_plants ? all_market_plant_keys :
                obj in mrkts ? all_market_keys : 
                obj in plants ? all_plant_keys : 
                all_keys : 
            arc_keys
        for obj in env_obj
    )
    #add keys to all nodes that place requests (non-source nodes = have a predecessor node)
    for node in nodes
        if node in nonsources
            push!(env_keys[node], :supplier_priority, :partial_fulfillment, :early_fulfillment)
        end
        if node in nonsinks
            push!(env_keys[node], :customer_priority)
        end
    end

    return env_keys
end

"""
    print_warnings(truncate_flag::Bool, roundoff_flag::Bool, replace_flag::Bool)

Print warning statements for modified stochastic parameters.
"""
function print_warnings(truncate_flag::Bool, roundoff_flag::Bool, replace_flag::Bool)
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
    check_bill_of_materials!(network::MetaDiGraph, n::Int)

Validate bill of material input at plant `n`.
"""
function check_bill_of_materials!(network::MetaDiGraph, n::Int)
    param = get_prop(network, n, :bill_of_materials)
    #create a sparse NamedArray from a BOM Dictionary
    if param isa Dict 
        i = unique(first.(keys(param))) #input materials
        j = unique(last.(keys(param))) #output materials
        # mats = union(keys(param)...) #all materials
        bom = NamedArray( #initialize named array
            spzeros(length(i),length(j)),
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
    @assert param isa NamedArray
    @assert (names(param,1) ⊆ mats) && (names(param,2) ⊆ mats) "Bill of material at node $n contains material names that have not been specified in the network metadata."
    @assert param isa NamedArray "The bill of materials at node $n must be a NamedArray."
end

"""
    check_make_to_order!(network::MetaDiGraph, n::Int)

Check list of make to order materials at plant `n`.
"""
function check_make_to_order!(network::MetaDiGraph, n::Int)
    if !(:make_to_order in keys(props(network,n)))
        set_prop!(network, n, :make_to_order, []) #default is no make to order materials
    else
        mto = get_prop(network, n, :make_to_order)
        @assert mto isa Vector "Make-to-order materials for node $n must be a `Vector`."
        mats = get_prop(network, :materials)
        remove = setdiff(mto, mats)
        @assert isempty(remove) "Some make-to-order materials at node $n are not listed in the system materials:\n$(join(remove,"\n"))."
        not_mto = filter(i -> !isproduced(network, n, i), mto)
        @assert isempty(not_mto) "Some make-to-order materials are not produced at node $n:\n$(join(not_mto,"\n"))."        
    end
end

"""
    set_default!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Material)

Set parameter defaults.
"""
function set_default!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Material)
    param_dict = get_prop(network, obj..., key)
    if key in [:inventory_capacity, :production_capacity] #default is uncapacitated
        set_prop!(network, obj, key, merge(param_dict, Dict(mat => Inf)))
    elseif key == :demand_frequency #demand at every period
        merge!(param_dict, Dict(mat => 1))
    elseif key == :supplier_priority #random ordering of supplier priority
        merge!(param_dict, Dict(mat => inneighbors(network, obj) |> pred -> obj in pred ? [obj] ∪ pred : pred)) #if producer, set itself as first supplier priority
    elseif key == :customer_priority #random ordering of customer priority
        merge!(param_dict, Dict(mat => sort(setdiff(outneighbors(network, obj),[obj]))))
    elseif key == :partial_fulfillment #allow partial fulfillment by default
        merge!(param_dict, Dict(mat => true))
    elseif key == :market_partial_fulfillment #allow partial fulfillment by default
        merge!(param_dict, Dict(mat => true))
    elseif key == :early_fulfillment #allow early fulfillment by default (before service lead time expires)
        merge!(param_dict, Dict(mat => true))
    elseif key == :market_early_fulfillment #allow early fulfillment by default (before service lead time expires)
        merge!(param_dict, Dict(mat => true))
    else #all others default to 0
        set_prop!(network, obj..., key, merge(param_dict, Dict(mat => 0.)))
    end
    param_dict = get_prop(network, obj..., key)

    return param_dict
end

"""
    check_supplier_priority(network::MetaDiGraph, obj::Int, mat::Material)

Validate supplier priority input.
"""
function check_supplier_priority(network::MetaDiGraph, obj::Int, mat::Material)
    param_dict = get_prop(network, obj, :supplier_priority)
    param = param_dict[mat]
    #if singleton is provided, put it in a Vector
    if param isa Number 
        set_prop!(network, obj, :supplier_priority, merge(param_dict, Dict(mat => [param])))
    end
    #check suppliers are connected to that node
    param_dict = get_prop(network, obj, :supplier_priority)
    unique!(param_dict[mat]) #remove any duplicates
    @assert param_dict[mat] ⊆ inneighbors(network, obj) "One or more suppliers listed in the supplier priority of node $obj is not a predecessor of that node."

    return param_dict
end
"""
    check_customer_priority(network::MetaDiGraph, obj::Int, mat::Material)

Validate customer priority input.
"""
function check_customer_priority(network::MetaDiGraph, obj::Int, mat::Material)
    param_dict = get_prop(network, obj, :customer_priority)
    unique!(param_dict[mat]) #remove any duplicates
    @assert Set(param_dict[mat]) == Set(setdiff(outneighbors(network, obj), [obj])) "The customer priority list for node $obj must contain every successor to that node."

    return param_dict
end


"""
    check_stochastic_variable!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Material)

Validate inputs for stochastic variables (lead time, demand quantity).
"""
function check_stochastic_variable!(network::MetaDiGraph, key::Symbol, obj::Union{Int, Tuple}, mat::Material)
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
    Parameter $key for material $mat at $obj must be a Sampleable distribution or a sigleton (non-negative) Vector.
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
    roundoff_flag = (key in [:lead_time,:service_lead_time]) && (param isa Distribution{T, Continuous} where T || (param isa Vector && !iszero(mod.(param,1))))
    #update param dict
    param_dict = get_prop(network, obj..., key) 

    return param_dict, replace_flag, truncate_flag, roundoff_flag
end

"""
    update_stochastic_parameter!(env::SupplyChainEnv, key::Symbol, obj::Union{Int, Tuple}, mat::Material, value::Any)

Allows updating the value of one of the stochastic parameters (key options: `:demand_distribution`, `:lead_time`, `:service_lead_time`)
along `obj` (arc or node) for material `mat`.
"""
function update_stochastic_parameter!(env::SupplyChainEnv, key::Symbol, obj::Union{Int, Tuple}, mat::Material, value::Any)
    network = env.network #get network
    truncate_flag, roundoff_flag, replace_flag = false, false, false #initialize flags
    param_dict = get_prop(network, obj..., key) #get parameter dictionary
    set_prop!(network, obj..., key, merge(param_dict, Dict(mat => value))) #update the value
    param_dict, replace_flag, truncate_flag, roundoff_flag = check_stochastic_variable!(network, key, obj, mat) #validate value
    print_warnings(truncate_flag, roundoff_flag, replace_flag) #print any warnings

    return param_dict
end

"""
    store_node_materials!(network::MetaDiGraph, plants::Vector)
Record which materials have non-zero capacity at each node (can be stored at that node).
This improves performance for simulations with many materials. Sort in materials in topological order 
    for plant nodes to account for any make-to-order intermediates.
"""
function store_node_materials!(network::MetaDiGraph, plants::Vector)
    materials = get_prop(network, :materials)
    for n in vertices(network)
        n_cap = filter(#get nonzero inventory capacities
            !iszero ∘ last,  
            get_prop(network, n, :inventory_capacity)
        )
        n_mats = collect(keys(n_cap)) ∩ materials
        if n in plants
            bom = get_prop(network, n, :bill_of_materials)
            mat_graph = material_graph(bom)
            mat_nodes = sort_topological(mat_graph)
            mat_node_names = [mat_graph[i,:name] for i in mat_nodes]
            n_mats = mat_node_names ∪ n_mats
        end
        set_prop!(network, n, :node_materials, n_mats)
    end
end 