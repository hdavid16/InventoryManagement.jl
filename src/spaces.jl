function action_space(env::SupplyChainEnv)
    num_products = length(env.materials)
    num_edges = ne(env.network)
    return [Interval{:closed,:open}(0,Inf) for i in num_products * num_edges]
end

function state(env::SupplyChainEnv)

end

function state_space(env::SupplyChainEnv)
    
end
