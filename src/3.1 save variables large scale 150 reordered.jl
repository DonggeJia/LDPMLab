using Distributed
addprocs(4)
@everywhere begin
    using Pkg
    Pkg.activate(".")  # Activate the current environment
    Pkg.add("JLD2")  # Ensure the package is added
end
@everywhere using JLD2, FileIO, Distances, Colors, GeometryBasics, LinearAlgebra, DelimitedFiles, PlotlyJS, PhysicalMeshes
@everywhere @load "variables for geometry-large50_reordered.jld2"
lolk = PlotlyJS.plot(scatter(x=[1,2], y=[2,3], mode="lines"))

rm("variables for geometry-large150_reordered.jld2")
using JLD2, SparseArrays, FileIO, Symbolics, SparseDiffTools, Multibreak, BenchmarkTools, Distances, Colors, GeometryBasics, LinearAlgebra, DelimitedFiles, PhysicalMeshes, StaticArrays, PlotlyJS
@load "variables for geometry-kupher2_less points2.jld2"
@load "variables for geometry-large50_reordered.jld2"
"C:/Users/DOJ14/OneDrive - University of Pittsburgh/PITT/Code/Julia/test_Mechanical/variables for geometry-large50_reordered.jld2"

using PlotlyJS
using Pkg
Pkg.update()
Pkg.instantiate()
Pkg.build("PlotlyJS")
Pkg.gc()
using Pkg
using Pkg
Pkg.Registry.add(RegistrySpec(name="General", url="https://mirrors.bfsu.edu.cn/julia/static/registries/General.git"))

Pkg.Registry.add(RegistrySpec(name="General", url="https://mirrors.bfsu.edu.cn/julia/static/registries/General.git"))
Pkg.add("JuliaZH")
ENV["HTTP_PROXY"] = "http://localhost:7890"
ENV["HTTPS_PROXY"] = "http://localhost:7890"

using JuliaZH
JuliaZH.generate_startup("BFSU")
@load "C:/Users/DOJ14/OneDrive - University of Pittsburgh/PITT/Code/Julia/test_Mechanical/variables for geometry-large30.jld2"

@save "shape sequence.jld2" shapesequen_save
@save "Rsequence.jld2" Rsequen_save
@save "pointsfinal.jld2" pointsfinal