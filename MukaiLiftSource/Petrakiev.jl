cd("MukaiLiftSource")
using Pkg; Pkg.activate("MukaiLiftP6")
using MukaiLiftP6
using HomotopyContinuation
using Oscar
using LinearAlgebra


S_start, L_start = make_start();
parametrized_system, l_start, A_rand = make_poly_system(S_start, L_start);

S_target_random1 = randn(ComplexF64,21)
S_target_random2 = randn(ComplexF64,21)
S_target_random3 = randn(ComplexF64,21)

@time result1 = HomotopyContinuation.solve(parametrized_system, l_start; compile = false, start_parameters=S_start, target_parameters=S_target_random1)

@time result2 = HomotopyContinuation.solve(parametrized_system, solutions(result1); compile = false, start_parameters=S_start_random1, target_parameters=S_target_random2)

@time result3 = HomotopyContinuation.solve(parametrized_system, solutions(result2); compile = false, start_parameters=S_start_random2, target_parameters=S_target_random3)

@time result4 = HomotopyContinuation.solve(parametrized_system, solutions(result3); compile = false, start_parameters=S_start_random3, target_parameters=S_target_random1)

# compare result1 and result4