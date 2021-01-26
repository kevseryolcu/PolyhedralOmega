#using PolyhedralOmega, LinearAlgebra, Test
using LinearAlgebra, Test, PolyhedralOmega

include("../src/SymbolicCone.jl")

@testset "macmahon" begin
    @test PolyhedralOmega.macmahon([1 -1; -2 1], [0, 0]) != nothing
end
