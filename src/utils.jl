"""
    isproduced(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})

Check if material `mat` is produced in node `n`.
"""
function isproduced(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})
    !in(n, x.producers) && return false #n is not a plant
    bom = get_prop(x.network, n, :bill_of_materials)
    !in(mat, names(bom,2)) && return false #mat is not produced at this plant
    raws = filter(k -> k < 0, bom[:,mat]) #names of raw materials
    isempty(raws) && return false #mat is not produced at this plant (no raw materials are converted to mat)

    return true
end

"""
    isconsumed(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})

Check if material `mat` isa consumed in node `n`.
"""
function isconsumed(x::SupplyChainEnv, n::Int, mat::Union{Symbol,String})
    !in(n, x.producers) && return false #n is not a plant
    bom = get_prop(x.network, n, :bill_of_materials)
    !in(mat, names(bom,1)) && return false #mat is not consumed at this plant
    prods = filter(k -> k < 0, bom[mat,:]) #names of products made from mat
    isempty(prods) && return false #mat is not consumed at this plant (no products are made from mat)

    return true
end

"""
    get_capacity_and_supply(
        bom::NamedArray, rmat_names::Vector, cmat_names::Vector, 
        capacities::Dict, supply_grp::GroupedDataFrame
    )

Get available capacity and material supply at producer.
"""
function get_capacity_and_supply(
    n::Int, mat::Union{Symbol,String}, bom::NamedArray, 
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