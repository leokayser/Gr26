include("../Utilities.jl")


R,ϕ,vrs,M = plückercoordinates(2,6,QQ);

gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt");
gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt");

A = reshape(gr_start_param[2:length(gr_start_param)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
B = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6
#B_normal_form = B*inv(B[1:7,:])


numerical_plücker = poly_to_fp(ϕ) # Function that evaluates the polynomials in ϕ into the input
Z = numerical_plücker.(gr_start_system); # 14 points in P^14

Γ = reshape(hcat([B\z for z in Z]...),7,14)

Γ_norm = normalize_SApoints(Γ)
Γ_norm[:,8:14]*transpose(Γ_norm[:,8:14])






####################

G = QQMatrix_to_ComplexF64(randomSApoints(14))

M1 = QQMatrix_to_ComplexF64(reshape(rand(49),7,7))
D = QQMatrix_to_ComplexF64(diagm(randn(14)))

Γ = M1*G*D
N = normalize_SApoints(Γ)[:,8:14]
transpose(N)*N 
Cayley_orth_to_skew([:,8:14])
S1 = Cayley_orth_to_skew(Γ_norm[:,8:14])

S1 + transpose(S1)


Γ = randomSApoints(14)
####################