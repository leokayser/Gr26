include("make_start_system.jl")
using Latexify

#using Random
#sk = skew_matrix(rand(-4:4,21)).entries
#Gamma = (rand(-1:1,7,7)*hcat( I + sk, I - sk)*Diagonal(rand([-1;1],14)))[:,randperm(14)]

Gamma = [-7 8 7 -6 -4 -4 -6 -2 8 -6 8 4 -7 -5; -7 -2 5 -3 9 -1 7 2 0 1 1 -4 5 3; -2 4 4 -5 -6 -3 -4 -4 6 -4 4 6 -3 -5; 0 -1 0 -7 -1 -7 -1 -5 -1 -3 5 5 3 5; 2 -9 -2 0 6 -2 6 -1 -7 8 -8 -1 8 6; 6 1 -8 6 1 8 3 5 -1 -7 5 -7 2 2; 1 -6 1 -2 3 -4 5 -4 -8 7 -5 4 7 7]
latexify(Gamma)

Lambda = self_duality_control(Gamma)
print([round(Lambda[i,i]/Lambda[1,1]) for i in 1:14])

Gamma_norm, A, Lambda_scaling = normalize_SApoints(Gamma)

norm(A*Gamma_norm - Gamma*Lambda_scaling, Inf)
P = Gamma_norm[:,8:14]
S = cayley_num(P)

TNF = [I+S I-S]

skew_to_vector(S)