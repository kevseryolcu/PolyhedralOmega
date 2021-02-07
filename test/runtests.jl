#using PolyhedralOmega, LinearAlgebra, Test
using LinearAlgebra, Test, PolyhedralOmega

@testset "macmahon" begin
    #@test PolyhedralOmega.macmahon([1 -1; -2 1], [0, 0]) != nothing
    include("../src/SymbolicCone.jl")
end
