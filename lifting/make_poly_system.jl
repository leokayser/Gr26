#include("../Utilities.jl")
include("make_start_system.jl")
include("union_system.jl")





S_start, φ_start = make_start();


@var a[1:21]
@var s[1:21]
@var x[1:7]
@var p[1:15]

S, s_oscar = PolynomialRing(QQ, ["s$i" for i=1:21])

# s -> skew matrix
#s_t = t*S_start + (1-t)*s_target

# skew matrix -> O (poor man's cayley using adjunct)
OO_t = fake_Cayley_polynomial(7,s_oscar);

# O -> Γ = [I_7 | O]
Γ_oscar = hcat(Matrix(identity_matrix(S,7)), OO_t );
println("oscar_to_HC_Q:")
Γ = []
Γ_sys = []
for i in 1:(length(Γ_oscar))
    Base.vcat!(Γ, Γ, oscar_to_HC_Q(Γ_oscar[i], s))
    println("$i/14 done.")

    Base.vcat!(Γ_sys, Γ_sys, oscar_to_HC_Q(Γ_oscar[i], s))
end
Γ_sys = [System(vcat(Γ[:,i],a), variables=[s;a]) for i=1:14]

# a -> φ_a (using φ_start[1:12,:] from before)

φ_a = vcat(φ_start[1:12,:], reshape(a,3,7))
φ_sys = System(φ_a*x, variables=[x;a])


#Z_sys = [HomotopyContinuation.compose(φ_sys,Γ_sys[i]) for i=1:14];

#Z = φ_a * Γ;

#Z_sys = System(reshape(Z,210), variables=[s[:];a[:]]);

# p_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], p) for i=1:15], variables=p)

plück_φ_sys = HomotopyContinuation.compose(plück_sys,φ_sys)  
total_sys = [ HomotopyContinuation.compose(plück_φ_sys, Γ_sys[i]) for i=1:14  ]
# total_sys is a vector of Systems. 
# We need to concatenate all the equations in a unique system keeping the same 'a' variables and different 'x'
F = union_systems(total_sys)



# f = p_j(Z_i) for i=1..14, j=1..15
#total_sys = [HomotopyContinuation.compose(plück_sys, Z_sys[i]) for i=1:14];
#length(system) #210



C = System(system_hc, variables = a[1:21], parameters = s[1:21])


##############