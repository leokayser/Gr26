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

A_rand = [randn(ComplexF64,15,7) for _ in eachindex(l) ]

L_l = L_start + sum(l[i]*A_rand[i] for i in eachindex(l))
L_sys = System(L_l*x, variables=[x;l])


# q_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

Q = plück_sys(expressions(L_sys));

Q_sys = System(Q, variables=[x;l]);
equations = vcat([Q_sys([Γ[:,i];l]) for i in 1:14]...);


parameterized_system = System(equations, variables=l, parameters=s);


l_start = zeros(ComplexF64,21)
#parameterized_system(l_start,S_start)   # Also sanity check

R = RandomizedSystem(parameterized_system, 21);

S_target = rand(ComplexF64,21)

println("Starting solve")

@time result = HomotopyContinuation.solve(R, l_start; start_parameters=S_start, target_parameters=S_target)

result.path_results

norm(parameterized_system(solutions(result)[1], S_target), Inf)
