# InventoryManagement.jl:

*Discrete-time simulation environment for Inventory Management in Supply Chain Networks.*

![](assets/logo.png)

## Overview

*InventoryManagement.jl* allows modeling a multi-period multi-product supply chain network under stochastic stationary demand and stochastic lead times. A supply network can be constructed using the following types of nodes:
- `Producers`: Nodes where inventory transformation takes place (e.g., raw materials are converted to intermediates or finished goods). Material transformation, including reactive systems with co-products, are modeled using [Bills of Materials](https://en.wikipedia.org/wiki/Bill_of_materials) (see [Model Inputs section](#node-specific-parameters)).
- `Distributors`: Nodes where inventory is stored and distributed (e.g., distribution centers).
- `Markets`: Nodes where end-customers place final product orders (e.g., retailer). 

These types of nodes can be used to model the following components of a supply network:
- `Manufacturing Plants`: Plants are modeled using at least two nodes joined by a directed arc:
  - Raw material storage node: Upstream `producer` node that stores raw materials and the plant's `bill of materials`.
  - Product storage node: Downstream `distributor` node that stores the materials produced at the plant.
  - Arc: the time elapsed between the consumption of raw materials and the production of goods (production time) is modeled with the arc lead time.
- `Distribution Centers`: DCs are modeled using `distributor` nodes.

Note: Any node can be marked as a `market` node to indicate that there is external demand for one or more materials stored at that node. This allows external demand at distribution centers or at manufacturing plants (either for raw materials or products).

The simplest network that can be modeled is one with a single retailer (`distributor`) with external demand (`market`) that is supplied by a warehouse (`distributor`). However, more complex systems can be modelled as well.

When defining a supply network, a `SupplyChainEnv` object is created based on system inputs and network structure. This object can then be used to simulate the inventory dynamics under a stochastic environment and a specified inventory management policy. Stochasticity is modeled in the external demand quantities at the `market` nodes and in the lead times between connected nodes in the network. However, deterministic values can also be used for external demand or lead times if desired. In each period of the simulation, a decision-maker can specify inventory replenishnment orders throughout the network (refered to as `actions`), which consist of the inventory quantities requested by each node to each immediate supplier for each material in the system. If no action is taken during the simulation, the inventory levels will eventually be depleted by the external demand. Depending on the system configuration, unfulfilled external or internal demand can be either backlogged or considered a lost sale.

The `SupplyChainEnv` can also potentially be used in conjunction with [ReinforcementLearning.jl](https://github.com/JuliaReinforcementLearning/ReinforcementLearning.jl) to train a Reinforcement Learning `agent` that places replenishment orders (`actions`) throughout the network.

This package generalizes and extends and the inventory management environment available in [OR-Gym](https://github.com/hubbs5/or-gym).

## Dependencies

*InventoryManagement.jl* relies primarily on the following packages:
- [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl): Tabulate results.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl): Define probability distributions for the lead times on the network arcs and the demand quantities at the market nodes.
- [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl): Define supply network topology
- [MetaGraphs.jl](https://github.com/JuliaGraphs/MetaGraphs.jl): Specify system parameters for each node or arc, or for the system as a whole.

## Installation

The package can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```julia
pkg> add InventoryManagement
```

For the master branch, run:
```julia
pkg> add https://github.com/hdavid16/InventoryManagement.jl
```

## Contact

**Author**: Hector D. Perez\
**Position**: Ph. D. Candidate @ Carnegie Mellon University\
**Email**: hdperez@cmu.edu\
**Year**: 2021
