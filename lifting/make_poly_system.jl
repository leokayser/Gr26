#include("../Utilities.jl")
include("make_start_system.jl")

S_start, L_start = make_start();

@var l[1:69]
@var s[1:21]
@var x[1:7]
@var q[1:15]

S, s_oscar = polynomial_ring(QQ, ["s$i" for i=1:21])

# s -> skew matrix

OscarSkewMat = skew_matrix(s_oscar)
Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]

Γ = [I+Smat I-Smat]

#A_rand = [randn(ComplexF64,15,7) for _ in eachindex(l) ];
#JLD.save("Random_matrices.jld", "data", A_rand )
using JLD
A_rand = JLD.load("Random_matrices.jld")["data"];  

L_l = L_start + sum(l[i]*A_rand[i] for i in eachindex(l));

L_sys = System(L_l*x, variables=[x;l]);


# q_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

Q = plück_sys(expressions(L_sys));

Q_sys = System(Q, variables=[x;l]);
equations = vcat([Q_sys([Γ[:,i];l]) for i in 1:14]...);


parametrized_system = System(equations, variables=l, parameters=s);

l_start = zeros(ComplexF64,length(l))
parametrized_system(l_start,S_start)   # Also sanity check

#@time R = RandomizedSystem(parametrized_system, 69);
#18473.318498 seconds (45.75 M allocations: 36.230 GiB, 0.23% gc time, 0.01% compilation time)

#=
S_target = randn(ComplexF64,21)


tnf_homotopy = parameter_homotopy(InterpretedSystem(parametrized_system);
    start_parameters = S_start, target_parameters = S_target);
tracker = Tracker(tnf_homotopy; options = TrackerOptions(parameters = FAST_TRACKER_PARAMETERS) )
result = track.(tracker, [l_start])
=#

println("Starting solve")

#@time result = HomotopyContinuation.solve(R, l_start; start_parameters=S_start, target_parameters=S_target)
#=

@time result = HomotopyContinuation.solve(
    parametrized_system, l_start; start_parameters=S_start,
    target_parameters=S_target, show_progress = true)

4771.640128 seconds (85.71 M allocations: 38.527 GiB, 0.08% gc time, 0.59% compilation time: <1% of which was recompilation)
Result with 1 solution
======================
• 1 path tracked
• 1 non-singular solution (0 real)
• random_seesd: 0xeee98a5a

result.path_results

1-element Vector{PathResult}:
 PathResult:
 • return_code → :success
 • solution → ComplexF64[-0.0592647754972386 + 0.7751087101315427im, -1.7725139883174452 + 0.23733801775977859im, -0.6460917178956368 + 0.10523009851410838im, -0.9342263367081752 + 0.6313277892169233im, 0.4991449361209814 + 0.5223585099896476im, -0.4030060846550946 + 1.0374715506460088im, 1.081979675428538 + 0.5907240719502352im, 0.5278479103606262 - 0.4934902432254344im, -0.25943321062084385 - 0.8441368086234682im, 0.6755946779669075 + 1.120391351798026im  …  0.1420331587427683 - 0.3215220677270636im, -0.8928647863133135 - 0.147063960980485im, 1.5657473829774033 - 0.5284368395846664im, -0.0031015347784415794 - 0.36794536804961264im, -1.1787553750331534 + 0.05393862147165947im, 0.5848935991138667 + 2.1717731026111977im, -0.701832739037763 + 0.5303619182739485im, -0.78366687719798 - 0.40826844654642896im, 1.0682398376061473 - 0.8145841536433973im, -0.4798529181072345 + 0.36161477022425537im]
 • accuracy → 9.2352e-17
 • residual → 4.0194e-13
 • condition_jacobian → 201.79
 • steps → 1412 / 0
 • extended_precision → false
 • path_number → 1


sol = solution(result[1])

norm(parametrized_system(sol, S_target), Inf)
 =#
