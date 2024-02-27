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
#OO_t = fake_Cayley_polynomial(7,s_oscar);
OscarSkewMat = skew_matrix(s_oscar)
Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]

Γ = [I+Smat I-Smat]



# O -> Γ = [I_7 | O]
Γ_oscar = hcat(Matrix(identity_matrix(S,7)), OO_t );
println("oscar_to_HC_Q:")
#Γ = []
Γ_sys = []
for i in 1:14
    #hc_polys = [oscar_to_HC_Q(f, s) for f in Γ_oscar[:,i]]
    #push!(Γ, hc_polys)
    push!(Γ_sys, System([Γ[:,i];a], variables=a, parameters=s))
    println("$i/14 done.")
end

# a -> φ_a (using φ_start[1:12,:] from before)

φ_a = vcat(φ_start[1:12,:], reshape(a,3,7))
φ_sys = System(φ_a*x, variables=[x;a])


#Z_sys = [HomotopyContinuation.compose(φ_sys,Γ_sys[i]) for i=1:14];

#Z = φ_a * Γ;

#Z_sys = System(reshape(Z,210), variables=[s[:];a[:]]);

# p_1..15 = plücker relations
plück_oscar = gens( grassmann_pluecker_ideal(2,6))

plück_sys = System([oscar_to_HC_Q(plück_oscar[i], p) for i=1:15], variables=p)

Q = plück_sys(expressions(φ_sys))

Q_sys = System(Q, variables=[x;a])
eqs = vcat([Q_sys([Γ[:,i];a]) for i in 1:14]...)

final_sys = System(eqs, variables=a, parameters=s)

plück_φ_sys = plück_sys ∘ φ_sys
sys_vector = Vector{AbstractSystem}()
for i in eachindex(Γ_sys)
    comp = plück_φ_sys ∘ Γ_sys[i]
    push!(sys_vector, comp);
    println("concat $i/14 done")
end
# total_sys is a vector of Systems. 
# We need to concatenate all the equations in a unique system keeping the same 'a' variables and different 'x'
F = union_system(sys_vector);

#φ_start = φ_start*inv(I+A)
a_start = reduce(vcat,φ_start[13:15,:])
φ_start == vcat(φ_start[1:12,:], reshape(a_start,3,7))
final_sys(a_start,S_start)

S_target = rand(ComplexF64,21)

@time result = HomotopyContinuation.solve(final_sys, a_start; start_parameters=S_start, target_parameters=S_target)
