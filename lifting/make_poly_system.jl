include("../Utilities.jl")
include("make_start_system.jl")

S_start, φ_start = make_start();


@var a[1:3,1:7]
@var s[1:21], s_target[1:21]
@var t

S, s_oscar = PolynomialRing(QQ, ["s$i" for i=1:21])

# s -> skew matrix
s_t = t*S_start + (1-t)*s_target

# skew matrix -> O (poor man's cayley using adjunct)
O_t = fake_Cayley_polynomial(7,s_oscar);

# O -> \Gamma = [I_7 | O]
Γ_oscar = hcat(Matrix(identity_matrix(S,7)), O_t )

# a -> \varphi_a (using \varphi_start[1:12,:] from before)
# Need to load \varphi_start from make_start_system

φ_a = vcat(φ_start[1:12,:], a)

# Z = \varphi_a * \Gamma
Γ =  [oscar_to_HC_Q(Γ_oscar_entry, s) for Γ_oscar_entry in Γ_oscar];

Z = φ_a * Γ;

# p_1..15 = plücker relations
Pluecker_rel = gens( grassmann_pluecker_ideal(2,6))
sys

# f = p_j(Z_i) for i=1..14, j=1..15
system = hcat([poly_to_fp(Pluecker_rel)(Z[:,i]) for i=1:14]...);
length(system) #210

##########
#C = System(system, variables = a, parameters = )



##############