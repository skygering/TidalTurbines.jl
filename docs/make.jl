using TidalTurbines
using Documenter

DocMeta.setdocmeta!(TidalTurbines, :DocTestSetup, :(using TidalTurbines); recursive=true)

makedocs(;
    modules=[TidalTurbines],
    authors="Skylar Gering, Quinn Early, and Weixuan Li",
    sitename="TidalTurbines.jl",
    format=Documenter.HTML(;
        canonical="https://skygering.github.io/TidalTurbines.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/skygering/TidalTurbines.jl",
    devbranch="main",
)
