module BHFunctions

# from basis.jl 
export generate_basis, build_index_dictionary

# from hamiltonian.jl
export occupation, initial_state, creation_ann, crea_crea_ann_ann, crea_crea_ann_ann, matrixoccupation,
generate_matrices,
matrix_of_crea_crea_ann_ann,
matrixofcreacreationannann,
H_theta

using LinearAlgebra #Uso de algebra lineal
#using Arpack #Encontrar eigen-valores
using SparseArrays #Uso de sparse-arrays 
#using ExponentialUtilities #Uso de matrices sparse en exponencial en un producto con un vector   

include("hamiltonian.jl")
include("basis.jl")

end
