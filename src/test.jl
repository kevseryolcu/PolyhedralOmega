module test

include("PolyhedralOmega.jl")

C = macmahon([2 2 6; 0 5 6], [1,2])
println("elim:", elimLastCoordinate(C))
# println("prim_v: ", prim_v([-2, 2, 4]))
println("prim: ", prim([-2 2 4; 0 5 6]));
# println("C: ", C)
# println("Flip res: ", flip(C))

end
