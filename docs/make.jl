push!(LOAD_PATH,"../src/")
using InventoryManagement
using Documenter
makedocs(
         sitename = "InventoryManagement.jl",
         modules  = [InventoryManagement],
         pages=[
                "Home" => "index.md",
                "Sequence of Events" => "events.md",
                "Inventory Policies" => "policies.md",
                "Supply Chain Model" => "model.md",
                "Examples" => "examples.md",
                "API" => "api.md"
               ])
deploydocs(;
    repo="github.com/hdavid16/InventoryManagement.jl",
)
