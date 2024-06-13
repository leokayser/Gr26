
using Oscar
using LinearAlgebra
using HomotopyContinuation
using DelimitedFiles
using JLD
using MukaiLiftP6

S_start, L_start = make_start();


parametrized_system, l_start = make_poly_system(S_start, L_start);
######## Solving for a random configuration ########

# S_target = randn(ComplexF64,21)
# @time result = HomotopyContinuation.solve(R, l_start; start_parameters=S_start, target_parameters=S_target)

####################################################