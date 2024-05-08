#include("../Utilities.jl")
include("make_start_system.jl")


S_start, L_start = make_start();


@var l[1:69]
@var s[1:21]
@var x[1:7]
@var q[1:15]

S, s_oscar = polynomial_ring(QQ, ["s$i" for i=1:21])

# s -> skew matrix
#s_t = t*S_start + (1-t)*s_target

OscarSkewMat = skew_matrix(s_oscar)
Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]

Γ = [I+Smat I-Smat]



# L -> L_a (using L_start[1:12,:] from before)

A_rand = [randn(ComplexF64,15,7) for _ in eachindex(l) ];

L_l = L_start + sum(l[i]*A_rand[i] for i in eachindex(l));

#indices = vcat([[i,1] for i in 2:15],[[i,2] for i in 2:9])

#L_l = Matrix{Any}(copy(L_start));
#for i in eachindex(l)
#    L_l[indices[i]...] = l[i]
#end




L_sys = System(L_l*x, variables=[x;l]);


# q_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

Q = plück_sys(expressions(L_sys));

Q_sys = System(Q, variables=[x;l]);
equations = vcat([Q_sys([Γ[:,i];l]) for i in 1:14]...);


parametrized_system = System(equations, variables=l, parameters=s);

#l_start = [L_start[indices[i]...] for i in eachindex(l)]

l_start = zeros(ComplexF64,length(l))
parametrized_system(l_start,S_start)   # Also sanity check

#R = RandomizedSystem(parameterized_system, 21);

S_target = randn(ComplexF64,21)

#=
tnf_homotopy = parameter_homotopy(InterpretedSystem(parametrized_system);
    start_parameters = S_start, target_parameters = S_target);
tracker = Tracker(tnf_homotopy; options = TrackerOptions(parameters = FAST_TRACKER_PARAMETERS) )
result = track.(tracker, [l_start])
=#

println("Starting solve")

#@time result = HomotopyContinuation.solve(R, l_start; start_parameters=S_start, target_parameters=S_target)
@time result = HomotopyContinuation.solve(
    parametrized_system, l_start; start_parameters=S_start,
    target_parameters=S_target, compile = false, show_progress = true)

result.path_results
sol = solution(result[1])

norm(parametrized_system(rand(ComplexF64,21), S_target), Inf)
