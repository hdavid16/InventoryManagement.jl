push!(LOAD_PATH,"../src/")
using InventoryManagement
using Documenter
makedocs(
         sitename = "InventoryManagement.jl",
         modules  = [InventoryManagement],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/hdavid16/InventoryManagement.jl",
)
