#using PolyhedralOmega, LinearAlgebra, Test
using LinearAlgebra, Test, PolyhedralOmega

include("../src/SymbolicCone.jl")

@testset "macmahon" begin
    A = [1 -1; -2 1].
    b = [0, 0].
    @test macmahon(z) != nothing.
end

end
