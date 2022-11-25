"""
    isproduced(net::Union{MetaDiGraph,SupplyChainEnv}, n::Int, mat::Material)

Check if material `mat` is produced in node `n`.
"""
function isproduced(net::MetaDiGraph, n::Int, mat::Material)
    !(:bill_of_materials in keys(props(net,n))) && return false #n is not a plant
    bom = get_prop(net, n, :bill_of_materials)
    !(mat in names(bom,2)) && return false #mat is not produced at this plant
    raws = filter(<(0), bom[:,mat]) #names of raw materials
    isempty(raws) && return false #mat is not produced at this plant (no raw materials are converted to mat)

    return true
end
isproduced(env::SupplyChainEnv, n, mat) = isproduced(env.network,n,mat)

"""
    isconsumed(net::Union{MetaDiGraph,SupplyChainEnv}, n::Int, mat::Material)

Check if material `mat` is consumed in node `n`.
"""
function isconsumed(net::MetaDiGraph, n::Int, mat::Material)
    !(:bill_of_materials in keys(props(net,n))) && return false #n is not a plant
    bom = get_prop(net, n, :bill_of_materials)
    !(mat in names(bom,1)) && return false #mat is not consumed at this plant
    prods = filter(<(0), bom[mat,:]) #names of products made from mat
    isempty(prods) && return false #mat is not consumed at this plant (no products are made from mat)

    return true
end
isconsumed(env::SupplyChainEnv, n, mat) = isconsumed(env.network,n,mat)

"""
    ismto(x::SupplyChainEnv, n::Int, mat::Material)

Check if material `mat` is make-to-order in node `n`.
"""
function ismto(x::SupplyChainEnv, n::Int, mat::Material)
    !(:make_to_order in keys(props(x.network,n))) && return false #n is not a plant
    mto = get_prop(x.network, n, :make_to_order)
    return mat in mto
end    

"""
    get_capacity_and_supply(
        x::SupplyChainEnv, n::Int, 
        mat::Material, bom::NamedArray, 
        rmat_names::Vector, cmat_names::Vector, 
        capacities::Dict
    )

Get available capacity and material supply at producer.
"""
function get_capacity_and_supply(
    x::SupplyChainEnv, n::Int, 
    mat::Material, bom::NamedArray, 
    rmat_names::Vector, cmat_names::Vector, 
    capacities::Dict
)
    #commit production at plant
    capacity = [] #get production capacity
    mat_supply = [] #store max capacity based on raw material consumption for each raw material
    isempty(rmat_names) && return [0],[0] #if material is not produced at the node, then return zero capacity and zero supply
    push!(capacity, capacities[mat]) #production capacity for that material
    for rmat in rmat_names #check raw material supply
        sup_pp = x.tmp[n,rmat,:on_hand] #supply of material involved in BOM
        push!(mat_supply, - sup_pp / bom[rmat,mat]) #only account for raw materials that are in the BOM
    end 
    for cmat in cmat_names #add capacity constraint for any co-products (scaled by stoichiometry)
        push!(capacity, capacities[rmat] / bom[cmat,mat])
    end

    return vcat(capacity, mat_supply)
end

"""
    sort_topological(net::MetaDiGraph)

Break self-loops and perform topological sort.
"""
sort_topological(net::MetaDiGraph) = topological_sort_by_dfs(remove_self_loops(net))

"""
    remove_self_loops(net::MetaDiGraph)

Create a copy of `net` with no self-loops.
"""
function remove_self_loops(net::MetaDiGraph)
    A = adjacency_matrix(net) #get adjacency_matrix
    A[diagind(A)] .= 0 #remove diagonal (self-loops)
    return SimpleDiGraph(A) #create new graph without self-loops
end

"""
    material_graph(bill_of_materials::NamedArray)

Create material graph based on a bill of materials.
"""
function material_graph(bill_of_materials::NamedArray)
    rmats = names(bill_of_materials, 1) #get input materials
    pmats = names(bill_of_materials, 2) #get output materials
    mats = union(rmats, pmats) #get all materials
    rids, pids, _ = findnz(sparse(bill_of_materials)) #get non-zero locs
    #create graph
    g = MetaDiGraph(length(mats)) 
    set_indexing_prop!(g, :name) #save node (material) names
    for n in vertices(g)
        set_prop!(g,n,:name,mats[n])
    end
    for (r,p) in zip(rids,pids)
        add_edge!(g, g[rmats[r], :name], g[pmats[p], :name])
    end

    return g
end

"""
    get_lead_time(x::SupplyChainEnv, arc::Edge, mat::Material)
"""
function get_lead_time(x::SupplyChainEnv, arc::Edge, mat::Material)
    lt_data = get_prop(x.network, arc, :lead_time_data)[mat]
    if !isnothing(lt_data)
        lt = lt_data[x.period]
    else
        lt = rand(get_prop(x.network, arc, :lead_time)[mat])
    end

    return ceil(Int,lt)
end