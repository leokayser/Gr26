# pkg> activate MukaiLiftP6
using Oscar
using LinearAlgebra
using HomotopyContinuation
using DelimitedFiles
using JLD
using MukaiLiftP6
using Latexify

S_start, L_start = make_start();
parametrized_system, l_start = make_poly_system(S_start, L_start);

#Sk0 = [0 3 1 2 1 1 2; -3 0 -1 0 -2 -3 4; -1 1 0 0 -1 -4 -2; -2 0 0 0 -3 2 2; -1 2 1 3 0 1 4; -1 3 4 -2 -1 0 0; -2 -4 2 -2 -4 0 0]
Gamma = [7 -2 6 -1 -6 1 -9 7 0 6 1 8 -3 7; -1 2 -5 -2 0 -4 3 -3 -4 -3 4 -2 4 -1; 1 4 -1 -5 -3 6 8 -1 -8 -3 5 1 -6 -8; 3 -6 4 -3 -4 6 0 5 8 2 3 2 -8 0; 1 -2 1 0 -4 2 2 3 4 -1 2 2 -2 -2; 0 -6 -5 6 3 7 -3 2 8 -7 -6 -3 -5 5; -3 3 -4 1 4 3 2 -3 -6 -4 -3 -4 -1 -2]
latexify(Gamma)

Lambda = self_duality_control(Gamma)
norm( Gamma*Lambda*transpose(Gamma), Inf)

Gamma_ONF, A, Lambda_scale = normalize_SApoints(Gamma)
A

norm(A*Gamma_ONF - Gamma*Lambda_scale, Inf)

P = Matrix{Float64}(Gamma_ONF[:,8:14])
latexify([I round.(P, digits=4)])
latexify(round.(Float64.(A/im), digits=4))

S = cayley_num(P)

Gamma_TNF = [I+S I-S]
latexify(Gamma_TNF)

S_target = skew_to_vector(S)

@time result = HomotopyContinuation.solve(parametrized_system, l_start; compile = false, start_parameters=S_start, target_parameters=S_target)
result.path_results
sol = solutions(result)[1]

L_hat = L_start + sum(sol[i]*A_rand[i] for i in eachindex(l))

print(latexify(round.(L_hat, digits=3)))

# parameterized_system(sol,S_target)

L = L_hat*(I+S)*inv(A)
q = poly_to_fp(plück_oscar)
maximum([ norm(q(L*Gamma[:,i]) ,Inf) for i =1:14])