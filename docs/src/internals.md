# Index

```@index
```

# Environment

```@docs
InventoryManagement.seed!
InventoryManagement.create_logging_dfs
```

# Inventory Policy

```@docs
InventoryManagement.is_valid
InventoryManagement.check_policy_inputs!
InventoryManagement.get_inventory_state
InventoryManagement.get_expected_consumption
InventoryManagement.calculate_reorder
```

# Reorder Actions

```@docs
InventoryManagement.initialize_inventories!
InventoryManagement.place_orders!
InventoryManagement.exit_order!
InventoryManagement.calculate_backlog
InventoryManagement.create_order!
InventoryManagement.fulfill_from_stock!
InventoryManagement.fulfill_from_production!
InventoryManagement.isproduced
InventoryManagement.get_capacity_and_supply
InventoryManagement.make_shipment!
InventoryManagement.log_unfulfilled_demand!
InventoryManagement.update_shipments!
InventoryManagement.enforce_inventory_limits!
InventoryManagement.update_inventories!
InventoryManagement.inventory_components
InventoryManagement.initialize_echelons!
InventoryManagement.update_echelons!
InventoryManagement.simulate_markets!
```

# Utilities
```@docs
InventoryManagement.show_action
InventoryManagement.check_inputs
```
