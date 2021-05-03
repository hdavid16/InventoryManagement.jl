function action_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_edges = ne(env.network)
    num_actions = num_products * num_edges
    ubound = []
    for n in vertices(env.network)
        set_prop!(env.network, n, :max_order, Dict(p => 0. for p in env.materials))
    end
    srcs = [n for n in nodes if isempty(inneighbors(network, n))]
    for (source, sink) in Iterators.material(srcs, env.markets)
        paths = yen_k_shortest_paths(env.network, source, sink, weights(env.network), typemax(Int)).paths
        #NOTE: TODO revise logic here if producers can be intermediate nodes!!!
        for path in paths
            capacity = get_prop(env.network, path[1], :production_capacity)
            for i in 2:length(path)
                top_max_order = get_prop(env.network, path[i-1], :max_order)
                if i > 2
                    top_init_inv = get_prop(env.network, path[i-1], :initial_inventory)
                end
                max_order = get_prop(env.network, path[i], :max_order)
                for p in env.materials
                    if i == 2
                        max_order[p] += capacity[p]
                    elseif i == 3
                        max_order[p] += top_init_inv[p] + top_max_order[p] * env.num_periods
                    else
                        max_order[p] += top_init_inv[p] + top_max_order[p]
                    end
                end
                set_prop!(env.network, path[i], :max_order, max_order)
            end
        end
    end
    for e in edges(env.network)
        max_order = get_prop(env.network, e.dst, :max_order)
        for p in env.materials
            push!(ubound, max_order[p])
        end
    end

    ClosedInterval.(zeros(num_actions), ubound)
end
