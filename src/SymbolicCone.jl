mutable struct SymbolicCone
    V::Array{Int64,2}
    q::Array{Float64,1}
    o::Array{Bool,1}
    sign::Int
    fp:: Array{Array{Int64},1}
    ratfun:: Any
    SymbolicCone(V::Array{Int64,2}, q::Array{Float64,1}, o::Array{Bool,1}, sign::Int=1) =
        new(V, q, o, sign)
end

function PrintSymbolicCone(C::SymbolicCone)
    println("\nSymbolic Cone:")
    println("V: \t", C.V)
    #display(C.V)
    println("q: \t", C.q)
    println("o: \t", C.o)
    println("sign:\t", C.sign)
    println("fp: \t", C.fp)
    println("Rational Function: ", C.ratfun, "\n\n")
end
