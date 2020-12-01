module MacMahonMethod
    using LinearAlgebra

    struct SymbolicCone
        V::Matrix{Int64}
        q::Vector{Int64}
        o::Vector{Bool}
        SymbolicCone(V::Matrix{T}, q::Vector{T}, o::Vector{Bool}) where T<:Int64 = new(V, q, o)
    end

    function macmahon(A::Matrix{Int64}, b::Vector{Int64})
        sizeA = size(A)
        Id = Matrix(1I, sizeA[2], sizeA[2])
        V = vcat(Id, A)
        q = append!(zeros(Int64, sizeA[2]), -b)
        o = zeros(Bool, sizeA[2])
        return SymbolicCone(V, q, o)
    end

    function flip(C::SymbolicCone)
        V = C.V
        s = 1
        for i = 1:size(V, 2)
            for j = 1:size(V, 1)
                if(V[j, i] != 0)
                    forward = (sign(V[j, i]) == 1)
                    println("forward: ", forward)
                    println("sign: ", sign(V[j, i]))
                    if(!forward)
                        s = s*(-1)
                        V[:, i] = (-1)*V[:, i]
                        C.o[i] = xor(C.o[i], true)
                    end
                    break
                end
            end
        end
        return (C, s)
    end

    function elimLastCoordinate(C::SymbolicCone)
        V = C.V
        n = size(V, 1)
        k = size(V, 2)
        q = C.q
        res = []
        if(q[n] < 0)
            II = filter(x -> V[n, x] > 0, Vector(1:k))
            println("Ay: ", II)
            Id = Matrix(1I, size(V, 1), size(V, 2))
            for j in II
                T = -1*Id*V[n, j]
                z = zeros(Int64, size(V, 1), size(V, 2))
                z[j, :] = V[n, :]
                T = T + z
                T[j, j] = -1
                println("T: ", T)
            end
        else
            II = filter(x -> V[n, x] < 0, Vector(1:k))
            println("Ay: ", II)
            Id = Matrix(1I, size(V, 1), size(V, 2))
            T = []
            for j in II
                T = -1*Id*V[n, j]
                z = zeros(Int64, size(V, 1), size(V, 2))
                z[j, :] = V[n, :]
                T = T + z
                T[j, j] = -1
                println("T: ", T)
            end
            println("TT: ", T)
            prim_res = prim(V)
            println("prim res: ", prim_res)
            C.V = T + prim_res
            println("T + prim_res: ", C.V)
        end

        return flip(C)

    end

    function prim_v(v::Array{Int64})
        d = abs(gcd(v))
        if(d == 1)
            return v
        else
            return transpose(reduce(hcat, map(vi -> floor(Int, vi/d), v))) #is this true?
        end
    end

    function prim(V::Matrix{Int64})
        return transpose(reduce(hcat, map(i -> prim_v(V[i, :]), Vector(1:size(V, 1)))))
    end

    C = macmahon([-1 -2 -3; 0 -5 -6], [1,2])
    println("C: ", C)
    println("Flip res: ", flip(C))
    elimLastCoordinate(C)
    println("prim_v: ", prim_v([-2, 2, 4]))
    println("prim: ", prim([-2 2 4; 0 5 6]));

end
