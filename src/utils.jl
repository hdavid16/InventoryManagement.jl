################################################################################
#                       DEPRECATED
################################################################################

"""
    material_conversion(net::MetaDiGraph, plants::Vector, mrkts::Vector)

Map the amount of raw materials required at each plant to produce one unit of finished good at a market.
This is used to update the echelon inventory position at the plants.

NOTE: There must be at most 1 path between each plant and market. Otherwise, there would be issues if the
bill of materials on one of the paths is different from the other paths for the same finished good.
"""
function material_conversion(net::MetaDiGraph, plants::Vector, mrkts::Vector)
    mats = get_prop(net, :materials)
    mat_conv = Dict(
        (plant,mrkt) => 
            path_material_conversion(net, plant, mrkt, plants)
        for (plant,mrkt) in Iterators.product(plants,mrkts)
    )

    return mat_conv
end

"""
    path_material_conversion(net::MetaDiGraph, plant::Int, mrkt::Int, plants::Vector)

Get the material conversion dictionary for all raw materials consumed at `plant` and products sold at `mrkt`.
"""
function path_material_conversion(net::MetaDiGraph, plant::Int, mrkt::Int, plants::Vector)
    #materials in network
    mats = get_prop(net, :materials)
    num_mats = length(mats) 
    plant_inv_cap = get_prop(net, plant, :inventory_capacity)
    mrkt_inv_cap = get_prop(net, mrkt, :inventory_capacity)
    plant_mats = keys(filter(i -> !iszero(i[2]), plant_inv_cap))
    mrkt_mats = keys(filter(i -> !iszero(i[2]), mrkt_inv_cap))
    #path between plant and market
    net_paths = yen_k_shortest_paths(net, plant, mrkt, weights(net), 100)
    isempty(net_paths.paths) && return Dict()
    #create network of materials
    mat_graph = MetaDiGraph(num_mats)
    set_indexing_prop!(mat_graph, :name)
    for (i,mat) in enumerate(mats)
        set_prop!(mat_graph, i, :name, mat)
    end
    for path in net_paths.paths, n in filter(i -> i in plants, path)
        bom = get_prop(net, n, :bill_of_materials)
        in_mats = names(bom,1)
        out_mats = names(bom,2)
        locs = findall(i -> i < 0, bom)
        for loc in locs
            in_mat = in_mats[loc[1]]
            out_mat = out_mats[loc[2]]
            stoich = bom[loc]
            src = mat_graph[in_mat,:name]
            dst = mat_graph[out_mat,:name]
            new_edge = add_edge!(mat_graph, src, dst, :weight, stoich)
            if !new_edge
                wt = get_prop(mat_graph, src, dst, :weight)
                if stoich < wt
                    set_prop!(mat_graph, src, dst, :weight, stoich)
                end
            end
            !new_edge && @warn """
            There is more than one way of producing $out_mat from $in_mat on the path between $plant and $mrkt. 
            The one that consumes the most $out_mat is used.
            """
        end
    end
    conversion_dict = Dict(
        (mraw,mprod) => 
            map_conversion(mat_graph, mraw, mprod)
        for (mraw,mprod) in zip(plant_mats, mrkt_mats)
    )

    return conversion_dict
end

"""
    map_conversion(mat_graph::MetaDiGraph, mraw::Symbol, mprod::Symbol)

Get the amount of `mraw` consumed to produce one unit of `mprod`.
"""
function map_conversion(mat_graph::MetaDiGraph, mraw::Symbol, mprod::Symbol)
    mat_paths = yen_k_shortest_paths(mat_graph, mat_graph[mraw,:name], mat_graph[mprod,:name], weights(mat_graph), 100)
    isempty(mat_paths.paths) && return 0
    length(mat_paths.paths) > 1 && @warn """
    There is more than one pathway for converting $mraw to $mprod. Only one of these is used.
    """
    path = mat_paths.paths[1]
    path_conv = [get_prop(mat_graph, path[k], path[k+1], :weight) for k in 1:length(path)-1] #stoichiometry for each material conversion on the path
    conv = -abs(prod(path_conv)) #multiply each stoichiometry and force to negative (consumption)
    return conv
end

"""
    material_conversion(net::MetaDiGraph)

Generate a material graph where the weights are the stoichiometry. 
Also generate a dictionary mapping the amount of each material consumed when one unit of product is procured.
NOTE: assumes that byproducts don't have their own primary pathways
      assumes that byproducts are only made in one of the pathways
      assumes a general bill of materials is passed for the whole network
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