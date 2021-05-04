# Examples

## Example 1

This example is for a 100 period simulation of a supply network with one plant (node 1) that supplies a retailer (node 3), with stochastic demand for product `:A`. Node 3, has an alternate supplier, which is a distribution center (node 2). Node 3 prefers replenishing from the plant, which has a lower lead time. A `(s,S)` reorder policy is used at the retailer. When the on-hand level for material `:A` is depleted at the plant, the plant begins transforming raw material `:B` into `:A`. There is limited raw material supply at the plant. When the raw material stocks-out, node 3 switches to node 2 for its supply.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex1.jl).*

![](assets/ex1_profit.png)
![](assets/ex1_position.png)

## Example 2

This example is for a 100 period simulation of a supply network with one warehouse (node 1) that supplies a retailer (node 2) with stochastic demand for product `:A`. A `(r,Q)` reorder policy is used at the retailer every 25 periods.

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex2.jl).*

![](assets/ex2_inventory.png)

## Example 3

This example shows how the simulator can be used for a chemical production system with co-production and material recycle. The system modeled is the batch plant described in [Kondili, et al. (1993)](https://www.sciencedirect.com/science/article/pii/009813549380015F?via%3Dihub), which can be modeled as the following supply network:

![](assets/ex3_schematic_drawio.png)

*See code [here](https://github.com/hdavid16/InventoryManagement.jl/blob/master/examples/ex3.jl).*

![](assets/ex3_intermediate_tanks.png)
![](assets/ex3_product_tanks.png)
![](assets/ex3_sales.png)
