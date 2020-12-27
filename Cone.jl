include("SymbolicCone.jl")

mutable struct Cone{T<:Number}
    S::SymbolicCone{T}
    sign::Int
    Cone{T}(S::SymbolicCone{T}, sign::Int) where {T<:Number} =
        new(S, sign)
end

function PrintCone(C::Cone)
    PrintSymbolicCone(C.S)
    println("sign: ", C.sign)
end
