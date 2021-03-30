include("PolyhedralOmega.jl")



#println("solve: ", solve([1 -1], [0]))
#println("solve: ", solve([1 -1; -2 1], [0, 0]))
println("solve: ", PolyhedralOmega.solve([-1 -1; 1 -1], [-4, 0]))
