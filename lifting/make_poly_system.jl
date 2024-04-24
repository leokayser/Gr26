#include("../Utilities.jl")
include("make_start_system.jl")


S_start, L_start = make_start();


@var l[1:21]
@var s[1:21]
@var x[1:7]
@var q[1:15]

S, s_oscar = polynomial_ring(QQ, ["s$i" for i=1:21])

# s -> skew matrix
#s_t = t*S_start + (1-t)*s_target

# skew matrix -> O (poor man's cayley using adjunct)
#OO_t = fake_Cayley_polynomial(7,s_oscar);
OscarSkewMat = skew_matrix(s_oscar)
Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]

Γ = [I+Smat I-Smat]



# L -> L_a (using L_start[1:12,:] from before)

L_l = vcat(L_start[1:12,:], reshape(l,3,7))
L_sys = System(L_l*x, variables=[x;l])


# q_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

Q = plück_sys(expressions(L_sys));

Q_sys = System(Q, variables=[x;l]);
equations = vcat([Q_sys([Γ[:,i];l]) for i in 1:14]...);

parameterized_system = System(equations, variables=l, parameters=s);


#L_start = L_start*inv(I+A)
l_start = reduce(vcat,L_start[13:15,:])
#L_start == vcat(L_start[1:12,:], reshape(l_start,3,7))  # Sanity check
# final_sys(l_start,S_start);                            # Also sanity check



#S_target = rand(ComplexF64,21)

#@time result = HomotopyContinuation.solve(parameterized_system, l_start; start_parameters=S_start, target_parameters=S_target, show_progress=false)
