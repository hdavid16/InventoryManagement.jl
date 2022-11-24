module InventoryManagement

using DataFrames, NamedArrays, SparseArrays
using Chain
using Random
using Distributions
using IntervalSets
using LinearAlgebra
using Graphs, MetaGraphs
import StatsBase: mean, std

const Material = Union{String, Symbol}

include("environment.jl")
include("actions.jl")
include("demand.jl")
include("financials.jl")
include("fulfillment.jl")
include("inventory.jl")
include("metrics.jl")
include("network.jl")
include("parameters.jl")
include("policy.jl")
include("spaces.jl")
include("utils.jl")
include("balance.jl")

export SupplyChainEnv, reset!, reward, is_terminated
export state, state_space, action_space, show_action, show_state
export reorder_policy, simulate_policy!
export connect_nodes!
export update_stochastic_parameter!
export isproduced, isconsumed, material_graph
export calculate_service_measures!
export inventory_balance, normal_base_stock
export MetaDiGraph, set_prop!, set_props!

end
