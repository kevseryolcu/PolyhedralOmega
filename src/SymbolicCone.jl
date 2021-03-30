mutable struct SymbolicCone
    V::Array{Int64,2}
    q::Array{Float64,1}
    o::Array{Bool,1}
    sign::Int
    fp:: Array{Array{Int64},1}
    ratfun:: Any
    SymbolicCone(V::Array{Int64,2}, q::Array{Float64,1}, o::Array{Bool,1}, fp::Array{Array{Int64},1}, ratfun::Any, sign::Int=1) =
        new(V, q, o, sign, fp, ratfun)
    SymbolicCone(V::Array{Int64,2}, q::Array{Float64,1}, o::Array{Bool,1}, sign::Int=1) =
        new(V, q, o, sign)
end


function Base.println(C::SymbolicCone)
    println("Symbolic Cone: \t", C.V, "\t", C.q, "\t", C.o, "\t", C.sign, "\n")
end

function Base.isequal()
    return true;
end
