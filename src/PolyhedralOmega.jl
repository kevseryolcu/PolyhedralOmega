module PolyhedralOmega
using LinearAlgebra
using IterTools
#using LazySet

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
    Iplus = filter(x -> V[n, x] > 0, Vector(1:k))
    Iminus = filter(x -> V[n, x] < 0, Vector(1:k))
    if (q[n] < 0 || (q[n] == 0 && size(Iplus, 1) <= size(Iminus, 1) && size(Iplus, 1) >= 1) )
        II = Iplus
    else
        res = [SymbolicCone(prim(V[1:n-1, :]), q[1:n-1], o[1:n-1])]
        II = Iminus
    end
    #Cprime if q[n] is zero

    for j in II
        L = []
        for i in 1:k
            if i == j
                if q[n] >= 0
                    push!(L, -1*V[:,j])
                 else
                    push!(L, V[:,j])
                 end
            else
                push!(L, (q[n] >= 0 ? 1 : -1)*( (V[n,i] * V[:,j]) + (-V[n,j]*V[:,i]) ))
            end
        end

        new_o = o[j] != false ? vcat(o[1:j-1], false, o[j+1:end]) : deepcopy(o)
        new_V = prim(hcat(L...))
        new_V = new_V[1:(size(new_V, 1)-1),:]

        new_q = map(+, q, (-q[n] / (V[n, j])) * V[:, j])
        new_q = new_q[1:n-1]
        new_cone = SymbolicCone(new_V, new_q, new_o)

        flip_res = flip(new_cone)
        push!(res, flip_res)
    end
    #######
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
    #println("elc res 1: ", res)
    res2 = []
    #println("K:", k)
    for i = 1:(k-1)#becuse we call eliminate last coordinate before
        for j = 1:length(res)
            eres = elimLastCoordinate(res[j])
            #println("eres: ", eres)
            append!(res2, eres)
            #println("res2: ", res2)
        end
        res = res2
        res2 = []
    end
    println("e res: ", res)
    return res
end

function enumerateFundamentalParallelePiped(C::SymbolicCone)
    SMFRes = SmithNormalForm.smith(C.V)
    S = SmithNormalForm.diagm(SMFRes)
    Uinv = inv(SMFRes.S)
    Winv = inv(SMFRes.T)
    dimension = size(C.V, 2) # num of rows
    ambientDimension = size(C.V, 1) # num of cols
    #println("C.V: ", C.V)
    #println("S: ", S, "\nUinv: ", Uinv, "\nWinv: ", Winv)
    #println("size 1: ", size(S,1), "\nsize 2: ", size(S,2))

    diagonals = Int64[]
    for i in 1:dimension
        if(i <= size(S,1) && i <= size(S,2))
            #println("i: ", i)
            push!(diagonals, S[i,i])
        end
    end
    #println("Diagonals: ", diagonals)

    lastDiagonal = diagonals[end]

    # sprime = [Integer(sk / si) for si in s]
    primeDiagonals = Int64[]
    for d in diagonals
        push!(primeDiagonals, Int64(lastDiagonal/d))
    end
    #println("Prime diagonals: ", primeDiagonals)

    # qhat = Uinv * q
    apex = C.q
    #println("q: ", apex, "\nV: ", C.V)

    qhat = Uinv*apex
    #println("qhat: ", qhat)

    # Wprime
    Wprime = [Winv[j,i]*primeDiagonals[i] for i = 1:dimension, j = 1:dimension] #Winv * primeDiagonals
    #println("wprime: ", Wprime)

    # qtrans
    qtrans = [sum([-Wprime[j,i] * qhat[i] for i = 1:dimension]) for j = 1:dimension]
    #println("qtrans: ", qtrans)

    #qfrac
    qfrac = [qtrans[i] - floor.(Int, qtrans[i]) for i = 1:dimension]
    #println("qfrac: ", qfrac)

    #qint
    qint = [ floor.(Int, qi) for qi in qtrans ]
    #println("qint: ", qint)

    #qsummand
    qsummand = [Int64(qi) for qi in (lastDiagonal*apex + C.V*qfrac) ]
    #println("qsummand", qsummand)

    #openness
    openness = [ (qfrac[j] == 0 ? C.o[j] : 0) for j in 1:dimension]
    #println("openness: ", openness)

    #bigP
    #res1 = [[1:1:diagonals[i];] for i= 1:dimension]
    #println("res1: ", res1)

    # CartesianProduct( *[xrange(s[i]) for i in 1:k] )
    L = []
    P = []
    for v in IterTools.product([1:diagonals[i] for i in 1:dimension]...)
        push!(P, v)
        innerRes = []
        j = 1
        for qj in qint
            inner = 0
            i = 1
            for vi in v
                inner += Wprime[i,j] * vi
                i += 1
            end
            inner += qj
            inner = inner % lastDiagonal

            if inner == 0 && C.o[j]
                inner = lastDiagonal
            end
            append!(innerRes, inner)
            j += 1
        end
        #println("inneerres: ", innerRes)

        outerRes = []
        for l in 1:ambientDimension
            outer = 0
            j = 1
            for innerResi in innerRes
                outer += C.V[l,j] * innerResi
                j += 1
            end
            append!(outerRes, outer) # outerRes is an integral vector
        end
        #println("outerRes: ", outerRes)
        push!(L, tuple(collect( ((ai + bi) / lastDiagonal) for (ai,bi) in collect(zip(outerRes, qsummand)) )))
    end
    return L
end


function solve(A::Matrix{Int64}, b::Vector{Int64})
    C = macmahon(A, b)
    ListOfSymbolicCones = eliminateCoordinates(C, size(C.V, 1) - size(C.V, 2))
    ListOfFundPP = []

    for cone in ListOfSymbolicCones
        PrintSymbolicCone(cone)
        append!(ListOfFundPP, enumerateFundamentalParallelePiped(cone))
    end

    #transformIntegral([1 -1])
    for f in ListOfFundPP
        for p in f
            println("p: ", collect(p))
        end
    end
    #println("FPP: ", ListOfFundPP)
    # list of rational_function

end


#C = macmahon([1 -1; -2 1], [0, 0])


#PrintSymbolicCone(C)

#println("C: ", C)
#println("Flip res: ", flip(C))
#res = elimLastCoordinate(C)
#println("res: ", res)
#PrintCone(res[1])

# println("eliminate cordinates result:\n", eliminateCoordinates(C, size(C.V, 1) - size(C.V, 2)))
# println("EFP: ", enumerateFundamentalParallelePiped(C))
#println("solve: ", solve([1 -1], [0]))
println("solve: ", solve([1 -1; -2 1], [0, 0]))
V = [1 1; 0 3]
q = [0.0, 0.0]
o = [false, false]
s = SymbolicCone(V, q, o, 1)
#PrintSymbolicCone(s)
#res = enumerateFundamentalParallelePiped(s)
#println("resfpp: ", res)
#println("prim_v: ", prim_v([-2, 2, 4]))
#println("prim: ", prim([-2 2 4; 0 3 6]))
end
