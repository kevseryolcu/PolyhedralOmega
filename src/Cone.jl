include("SymbolicCone.jl")

mutable struct Cone
    S::SymbolicCone
    sign::Int
    Cone(S::SymbolicCone, sign::Int) =
        new(S, sign)
end

function PrintCone(C::Cone)
    PrintSymbolicCone(C.S)
    println("sign: ", C.sign)
end
