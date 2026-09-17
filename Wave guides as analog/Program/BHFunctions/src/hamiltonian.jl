function occupation(i, state_v) 
     if state_v[i] >= 1  
        return state_v[i]  
    else 
        return 0
    end  
end

# Estado inicial
function initial_state(v, D, M)
    vector = zeros(D)
    counter_position = 1;
    for i in v 
        if(i[Int(M/2)]==1 && i[Int(M/2 + 1)]==1) 
            vector[counter_position] = 1
            return vector
        end 
        counter_position+=1
    end  
    return "Deberia checar los sitios pares"
end

function creation_ann(k,j,v)     
    if(v[j] >= 1) 
        b = copy(v)
        value = sqrt(((b[k])+1)*b[j])
        b[j] -= 1; 
        b[k] += 1;    
        return [value , b]   
    else 
        return 0; 
    end 
end   

function crea_creation_ann_ann(k,j,v)     
    if(v[j] >= 1 && v[k] >= 1) 
        b = copy(v)   
        value = sqrt(b[j])
        b[j] -= 1 
        if(b[k] >= 1)  
            value *= sqrt(b[k])
            b[k] -= 1; 
            value *= sqrt((b[k])+1)
            b[k] += 1;
            value *= sqrt((b[j])+1)
            b[j] += 1;
            return [value , b] 
        else 
            return 0; 
        end
    else 
        return 0; 
    end 
end   

function crea_crea_ann_ann(k,c,j,g,v)     
    if(v[j] >= 1 && v[g] >= 1) 
        b = copy(v) 
        value = sqrt(b[g]) 
        b[g] -= 1 
        if(b[j] >= 1) 
            value *= sqrt(b[j]) 
            b[j] -= 1 
            value *= sqrt((b[c])+1) 
            b[c] += 1 
            value *= sqrt((b[k])+1) 
            b[k] += 1
            return [value , b]  
        else 
            return 0; 
        end
    else 
        return 0; 
    end 
end      

# ---------------------------------------------------------
# Generación de matrices y diccionarios
# ---------------------------------------------------------

# 
function matrixoccupation(k, v, D)
    diag_elements = [occupation(k, v[j]) for j in 1:D] 
    return sparse(1:D, 1:D, diag_elements, D, D)  
end 

# 
function generate_matrices(v, D, M)
    matricesdeoc = Vector{SparseMatrixCSC{Complex{Float32}, Int32}}(undef, M) 
    for d in 1:M
        matricesdeoc[d] = matrixoccupation(d, v, D)
    end
    return matricesdeoc
end
#
function matrix_of_crea_crea_ann_ann(i, d, k, g, D, M, v, DictOfInds)
    mat = zeros(D,D)
    for j in 1:D 
        if crea_crea_ann_ann(i, d, k, g, v[j])== 0 
           continue      
        end
        mat[DictOfInds[tag(crea_crea_ann_ann(i, d, k, g, v[j])[2],M)], j] = crea_crea_ann_ann(i, d, k, g, v[j])[1] 
    end    
    return mat
end

# 
function matrixofcreacreationannann(i, k, D, M, v, DictOfInds)
    mat = zeros(D,D)
    for j in 1:D 
        if crea_creation_ann_ann(i,k,v[j]) == 0 
           continue      
        end
        mat[DictOfInds[tag(crea_creation_ann_ann( i, k , v[j] )[2],M)], j] = crea_creation_ann_ann(i, k , v[j])[1] 
    end    
    return mat
end

# ---------------------------------------------------------
# Hamiltoniano
# ---------------------------------------------------------

# REORDENAMOS ARGUMENTOS: Los que NO tienen valor por defecto van PRIMERO.
function H_theta(v, M, D, DictOfInds, Mocup, θ::Float64, U::Float64=0.0, g::Float64=0.0, J::Float64=1.0)
    
    # Interaction term  
    Hint = spzeros(D,D)
    for d in 1:M 
        Hint += Mocup[d]*(Mocup[d] - I) 
    end  
    Hint = (U/2)*Hint 
    
    # Cavity term
    Hcav = spzeros(D,D)  
    V = g/M 
    for s in 1:M, u in 1:M 
        Hcav += (-1)^(s-u)*Mocup[s]*Mocup[u] 
    end 
    Hcav = -V*Hcav  
    
    # Potential part 
    Htheta = spzeros(Complex{Float64}, D, D)  
    for j in 1:D, k in 1:M, m in 1:M    
        if creation_ann(k,m,v[j]) == 0 
            continue     
        end 
        pos = DictOfInds[tag(creation_ann( k, m , v[j] )[2], M)] 
        if ( abs(k-m)==1 )    
            if(k-m==1)  
                Htheta[j, pos] = exp(-1im*θ*occupation(m,creation_ann( k, m , v[j])[2])) 
                continue
            end 
            if(k-m==-1)  
                Htheta[j, pos] = exp(1im*θ*occupation(k,v[j])) 
                continue
            end 
        end 
    end      
    
    H = spzeros(Complex{Float32}, D, D);  
    Hkin = spzeros(Complex{Float32}, D, D);  
    
    for j in 1:D, k in 1:M, m in 1:M    
        if creation_ann(k,m,v[j]) == 0 
            continue     
        end
        if (abs(k-m)==1)             
            if(DictOfInds[tag(creation_ann( k, m , v[j] )[2],M)] > j) 
                Hkin[DictOfInds[tag( creation_ann( k, m , v[j])[2] ,M)], j] = -J*(creation_ann( k, m , v[j])[1])   
            end
        end      
    end      
    
    Hkin = Symmetric(Hkin, :L) 
    H = Hkin.*(Htheta) .+ Hint .+ Hcav
    
    return H
end
