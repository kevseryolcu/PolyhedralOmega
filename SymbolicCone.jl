mutable struct SymbolicCone{T<:Number}
    V::Matrix{T}
    q::Vector{Float64}
    o::Vector{Bool}
    SymbolicCone{T}(V::Matrix{T}, q::Vector{Float64}, o::Vector{Bool}) where {T<:Number} =
        new(V, q, o)
end

function PrintSymbolicCone(C::SymbolicCone)
    println("\nSymbolic Cone:")
    println("V:")
    display(C.V)
    println("q: ", C.q)
    println("o: ", C.o)
end
