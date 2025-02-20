#cd("MukaiLiftSource")
using Pkg; Pkg.activate("MukaiLiftP6")
using MukaiLiftP6
using HomotopyContinuation
using Oscar
using LinearAlgebra


S_start, L_start = make_start();
parametrized_system, l_start, A_rand = make_poly_system(S_start, L_start);

#S_target_random1 = randn(ComplexF64,21)
#S_target_random2 = randn(ComplexF64,21)
#S_target_random3 = randn(ComplexF64,21)

HomotopyContinuation.write_parameters("S_random1.txt", S_target_random1)
HomotopyContinuation.write_parameters("S_random2.txt", S_target_random2)
HomotopyContinuation.write_parameters("S_random3.txt", S_target_random3)

@time result_1a = HomotopyContinuation.solve(parametrized_system, l_start; compile = false, start_parameters=S_start, target_parameters=S_target_random1)
# 2726.222327 seconds
sol_1a = solutions(result_1a)
HomotopyContinuation.write_parameters("Solution_1a.txt", sol_1a[1])
#sol_1a = HomotopyContinuation.read_parameters("Solution_1a.txt")


@time result2 = HomotopyContinuation.solve(parametrized_system, sol_1a; compile = false, start_parameters=S_target_random1, target_parameters=S_target_random2)
#3873.805867 seconds
sol2 = solutions(result2)
HomotopyContinuation.write_parameters("Solution2.txt", sol2[1])
#sol_2 = HomotopyContinuation.read_parameters("Solution_2.txt")


@time result3 = HomotopyContinuation.solve(parametrized_system, sol2; compile = false, start_parameters=S_target_random2, target_parameters=S_target_random3)
#3109.211078 seconds 
sol3 = solutions(result3)
HomotopyContinuation.write_parameters("Solution3.txt", sol3[1])
#sol_3 = HomotopyContinuation.read_parameters("Solution_3.txt")


@time result_1b = HomotopyContinuation.solve(parametrized_system, sol3; compile = false, start_parameters=S_target_random3, target_parameters=S_target_random1)
#3677.871877 seconds
sol_1b = solutions(result_1b)
HomotopyContinuation.write_parameters("Solution_1b.txt", sol_1b[1])
#sol_1b = HomotopyContinuation.read_parameters("Solution_1b.txt")

sol_1a == sol_1b # compare this 

L_tilde_1a = L_start + sum(sol_1a[1][i]*A_rand[i] for i in eachindex(A_rand)) 
L_tilde_1b = L_start + sum(sol_1b[1][i]*A_rand[i] for i in eachindex(A_rand)) 

# L_tilde_1a == L_tilde_1b up to SO(7)???


# solve start system again
@time result_0 = HomotopyContinuation.solve(parametrized_system, sol2; compile = false, start_parameters=S_target_random2, target_parameters=S_start)
sol_0 = solutions(result_0)
HomotopyContinuation.write_parameters("Solution_0.txt", sol_0[1])
#sol_0 = HomotopyContinuation.read_parameters("Solution_0.txt")

L_tilde_0 = L_start + sum(sol_0[i]*A_rand[i] for i in eachindex(A_rand)) 


#####################################
########### Slicing Test ############
#####################################

S1 = vector_to_skew(S_target_random1)
Γ1 = vcat(I , cayley(S1) )
#L = L_tilde*(I+S)*inv(A)

## create system with pluecker relations
@var q[1:15]
plück_oscar = gens( grassmann_pluecker_ideal(2,6))
plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

maximum([ norm(plück_sys(L_tilde_1a * Γ1[:,i]) ,Inf) for i =1:14]) < 1e-10
maximum([ norm(plück_sys(L_tilde_1b * Γ1[:,i]) ,Inf) for i =1:14]) < 1e-10