mutable struct SymbolicCone
    V::Array{Int64,2}
    q::Array{Float64,1}
    o::Array{Bool,1}
    sign::Int
    SymbolicCone(V::Array{Int64,2}, q::Array{Float64,1}, o::Array{Bool,1}, sign::Int=1) =
        new(V, q, o, sign)
end

function PrintSymbolicCone(C::SymbolicCone)
    println("\nSymbolic Cone:")
    println("V:")
    display(C.V)
    println("q: ", C.q)
    println("o: ", C.o)
    println("sign: ", C.sign)
end
