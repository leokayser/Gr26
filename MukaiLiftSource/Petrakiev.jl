cd("MukaiLiftSource")
using Pkg; Pkg.activate("MukaiLiftP6")
using MukaiLiftP6
using HomotopyContinuation



S_start, L_start = make_start();
parametrized_system, l_start, A_rand = make_poly_system(S_start, L_start);

S_target_random = randn(ComplexF64,21)

@time result = HomotopyContinuation.solve(parametrized_system, l_start; compile = false, start_parameters=S_start, target_parameters=S_target_random)

