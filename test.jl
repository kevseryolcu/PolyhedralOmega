module test
    include("macmahonmethod.jl")

    C = macmahon([2 2 6; 0 5 6], [1,2])
    # println("C: ", C)
    # println("Flip res: ", flip(C))
    println("elim:", elimLastCoordinate(C))
    # println("prim_v: ", prim_v([-2, 2, 4]))
    println("prim: ", prim([-2 2 4; 0 5 6]));
end
