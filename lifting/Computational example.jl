#include("make_start_system.jl")
include("make_poly_system.jl")
using Latexify

#=
done = false
while !done
    sk = skew_matrix(rand(-4:4,21)).entries
    Gamma = rand(-1:1,7,7)*hcat( I + sk, I - sk)*Diagonal(rand([-1,1],14))
    done = true
    for a in Gamma
        if abs(a)>9
            done = false
        end
    end
end
=#

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

@time result = HomotopyContinuation.solve(parameterized_system, l_start; start_parameters=S_start, target_parameters=S_target)
#L_tilde = vcat(L_start[1:12,:], reshape(solutions(result)[1],3,7))
#print(latexify(round.(L_tilde, digits=3)))
sol = solutions(result)[1]
parameterized_system(sol,S_target)
parameterized_system(l_start,S_start)
result