#Generamos la base para N partículas y M sitios
function generate_basis(N, M)
    D = prod(max(N, M):N+M-1) ÷ prod(1:min(N, M))
    basis = [zeros(Int, M) for _ in 1:D]
    basis[1][1] = N
    for t = 2:D
        if basis[t-1][M] != 0
            k = M - 1
        else
            k = M
        end
        while k > 0 && basis[t-1][k] == 0
            k -= 1
        end
        @views basis[t][1:k-1] .= basis[t-1][1:k-1] 
        basis[t][k] = basis[t-1][k] - 1
        basis[t][k+1] = N - sum(@view(basis[t][1:k])) 
        end 
    return basis  
end  

# Función de números primos (Zhang)
p(i) = 100i + 3 

# Función tag 
function tag(b, M)
    tag_val = 0.0
    for i in 1:M
        tag_val += (√p(i)) * b[i] 	
    end
    return tag_val 
end

# Generamos el diccionario que mapea un estado a su índice en la matriz
function build_index_dictionary(basis, M)
    T = Float64[]
    for i in 1:length(basis) 
        push!(T, tag(basis[i], M)) 
    end 
    inds = sortperm(T)
    Tsorted = T[inds]
    return Dict(Tsorted .=> inds) 
end
