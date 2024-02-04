using LDPMLab
using Documenter

DocMeta.setdocmeta!(LDPMLab, :DocTestSetup, :(using LDPMLab); recursive=true)

makedocs(;
    modules=[LDPMLab],
    authors="Dongge Jia",
    sitename="LDPMLab.jl",
    format=Documenter.HTML(;
        canonical="https://DonggeJia.github.io/LDPMLab.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "Tutorials.md",
        "Example" => ["Mechanical response" => "mechanical_response.md",
            "Mass transport" => "mass_transport.md",
            "Bar-Reinforced LDPM" => "Embedded_bar_Reinforcement.md"],
        "Real applications" => "Real applications.md",
        "API" => "API.md",
        "References" => "references.md",
        "Contributing" => "contributing_guide.md"
    ])

deploydocs(;
    repo="github.com/DonggeJia/LDPMLab.jl",
    devbranch="main",
)
