module BHFunctions


using LinearAlgebra #Uso de algebra lineal
using Arpack #Encontrar eigen-valores
using SparseArrays #Uso de sparse-arrays 
using ExponentialUtilities #Uso de matrices sparse en exponencial en un producto con un vector   

include("hamiltonian.jl")
include("basis.jl")

export H_theta, generate_basis
end


