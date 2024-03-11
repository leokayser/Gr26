#include("../Utilities.jl")
include("make_start_system.jl")



S_start, φ_start = make_start();


@var a[1:21]
@var s[1:21]
@var x[1:7]
@var p[1:15]

S, s_oscar = PolynomialRing(QQ, ["s$i" for i=1:21])

# s -> skew matrix
#s_t = t*S_start + (1-t)*s_target

# skew matrix -> O (poor man's cayley using adjunct)
#OO_t = fake_Cayley_polynomial(7,s_oscar);
OscarSkewMat = skew_matrix(s_oscar)
Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]

Γ = [I+Smat I-Smat]



# a -> φ_a (using φ_start[1:12,:] from before)

φ_a = vcat(φ_start[1:12,:], reshape(a,3,7))
φ_sys = System(φ_a*x, variables=[x;a])


# p_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], p) for i=1:15], variables=p)

Q = plück_sys(expressions(φ_sys));

Q_sys = System(Q, variables=[x;a]);
eqs = vcat([Q_sys([Γ[:,i];a]) for i in 1:14]...);

final_sys = System(eqs, variables=a, parameters=s);


#φ_start = φ_start*inv(I+A)
a_start = reduce(vcat,φ_start[13:15,:])
#φ_start == vcat(φ_start[1:12,:], reshape(a_start,3,7))  # Sanity check
# final_sys(a_start,S_start);                            # Also sanity check

S_target = rand(ComplexF64,21)

@time result = HomotopyContinuation.solve(final_sys, a_start; start_parameters=S_start, target_parameters=S_target)
