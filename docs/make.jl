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
    ],
)

deploydocs(;
    repo="github.com/DonggeJia/LDPMLab.jl",
    devbranch="main",
)
