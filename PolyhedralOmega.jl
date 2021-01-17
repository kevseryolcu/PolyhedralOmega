module PolyhedralOmega
using LinearAlgebra

include("SymbolicCone.jl")
include("./SmithNormalForm/src/SmithNormalForm.jl")

function macmahon(A::Matrix{Int64}, b::Vector{Int64})
    sizeA = size(A)
    Id = Matrix(1I, sizeA[2], sizeA[2])
    V = vcat(Id, A)
    q = append!(zeros(Float64, sizeA[2]), -b)
    o = zeros(Bool, sizeA[2])
    return SymbolicCone(V, q, o)
end

function flip(C::SymbolicCone)
    V = C.V
    s = 1
    for i = 1:size(V, 2)
        for j = 1:size(V, 1)
            if (V[j, i] != 0)
                forward = (sign(V[j, i]) == 1)
                if (!forward)
                    s = s * (-1)
                    V[:, i] = (-1) * V[:, i]
                    C.o[i] = xor(C.o[i], true)
                end
                break
            end
        end
    end
    C.sign = s
    return C
end

function elimLastCoordinate(C::SymbolicCone)
    V = C.V
    n = size(V, 1)
    k = size(V, 2)
    q = C.q
    o = C.o
    res = []
    II = []
    # println("V: ", V)
    # println("q: ", q)
    # println("k: ", k)
    # println("vector: ", Vector(1:k))
    # println("filter1: ", filter(x -> V[n, x] > 0, Vector(1:k)))
    # println("filter2: ", filter(x -> V[n, x] <= 0, Vector(1:k)))
    #println("qn: ", q[n])

    if (q[n] < 0)
        append!(II, filter(x -> V[n, x] > 0, Vector(1:k)))
    else
        res = [SymbolicCone(prim(V[1:n-1, :]), q[1:n-1], o[1:n-1])]
        append!(II, filter(x -> V[n, x] < 0, Vector(1:k)))
    end
    #Cprime if q[n] is zero


    Id = Matrix(1I, size(V, 1), size(V, 2))
    T = []
    for j in II
        T = -1 * Id * V[n, j]
        z = zeros(Int64, size(V, 1), size(V, 2))
        z[j, :] = V[n, :]
        T = T + z
        T[j, j] = -1

        prim_res = prim(V)
        new_V = T + prim_res
        new_V = new_V[1:end.!=n, 1:end]

        new_o = o[j] != 0 ? o[:, j] + (0,) + o[j+1, :] : deepcopy(o)

        new_q = map(+, q, (-q[n] / (V[n, j])) * V[:, j])
        new_q = new_q[1:n-1]

        new_cone = SymbolicCone(new_V, new_q, new_o)

        flip_res = flip(new_cone)
        push!(res, flip_res)
    end
    return res
end

function prim_v(v::Array{Int64})
    d = abs(gcd(v))
    if (d == 1)
        return v
    else
        return reduce(vcat, map(vi -> floor(Int, vi / d), v))
    end
end

function prim(V::Matrix{Int64})
    return Matrix{Int64}(transpose(reduce(hcat, map(i -> prim_v(V[i, :]), Vector(1:size(V, 1))))))
end

function eliminateCoordinates(C::SymbolicCone, k::Int)
    res = elimLastCoordinate(C)
    res2 = []
    #println("K:", k)
    for i = 1:(k-1)#becuse we call eliminate last coordinate before
        for j = 1:length(res)
            eres = elimLastCoordinate(res[j])
            #println("eres: ", eres)
            push!(res2, eres)
            #println("res2: ", res2)
        end
        res = res2
        res2 = []
    end
    #println("e res: ", res)
    return res
end

function enumerateFundamentalParallelePipe
end

C = macmahon([1 -1], [0])


#PrintSymbolicCone(C)

println("C: ", C)
# println("Flip res: ", flip(C))
#res = elimLastCoordinate(C)
#println("res: ", res)
#PrintCone(res[1])

println("eliminate cordinates result:\n", eliminateCoordinates(C, size(C.V, 1) - size(C.V, 2)))
#println("prim_v: ", prim_v([-2, 2, 4]))
#println("prim: ", prim([-2 2 4; 0 3 6]))

end
