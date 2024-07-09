using LDPMLab
using Test

@testset "LDPMLab.jl" begin
    # Write your tests here.

    @test typeof(LDPM.geometry_parameters) == Vector{Float64}
end

