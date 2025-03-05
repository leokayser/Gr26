#cd("MukaiLiftSource")
using Pkg; Pkg.activate("MukaiLiftP6")
using MukaiLiftP6
using HomotopyContinuation
using Oscar
using LinearAlgebra


S_start, L_start = make_start();
parametrized_system, l_start, A_rand = make_poly_system(S_start, L_start);

#=
S_target_random1 = randn(ComplexF64,21)
S_target_random2 = randn(ComplexF64,21)
S_target_random3 = randn(ComplexF64,21)

HomotopyContinuation.write_parameters("S_random1.txt", S_target_random1)
HomotopyContinuation.write_parameters("S_random2.txt", S_target_random2)
HomotopyContinuation.write_parameters("S_random3.txt", S_target_random3)
=#

#=
@time result_1a = HomotopyContinuation.solve(parametrized_system, l_start; compile = false, start_parameters=S_start, target_parameters=S_target_random1)
# 2726.222327 seconds
sol_1a = solutions(result_1a)
HomotopyContinuation.write_parameters("Solution_1a.txt", sol_1a[1])
#sol_1a = HomotopyContinuation.read_parameters("Solution_1a.txt")


@time result2 = HomotopyContinuation.solve(parametrized_system, sol_1a; compile = false, start_parameters=S_target_random1, target_parameters=S_target_random2)
#3873.805867 seconds
sol2 = solutions(result2)
HomotopyContinuation.write_parameters("Solution_2.txt", sol2[1])
#sol_2 = HomotopyContinuation.read_parameters("Solution_2.txt")


@time result3 = HomotopyContinuation.solve(parametrized_system, sol2; compile = false, start_parameters=S_target_random2, target_parameters=S_target_random3)
#3109.211078 seconds 
sol3 = solutions(result3)
HomotopyContinuation.write_parameters("Solution_3.txt", sol3[1])
#sol_3 = HomotopyContinuation.read_parameters("Solution_3.txt")


@time result_1b = HomotopyContinuation.solve(parametrized_system, sol3; compile = false, start_parameters=S_target_random3, target_parameters=S_target_random1)
#3677.871877 seconds
sol_1b = solutions(result_1b)
HomotopyContinuation.write_parameters("Solution_1b.txt", sol_1b[1])
#sol_1b = HomotopyContinuation.read_parameters("Solution_1b.txt")


# solve start system again
@time result_0 = HomotopyContinuation.solve(parametrized_system, sol2; compile = false, start_parameters=S_target_random2, target_parameters=S_start)
#3790.035672 seconds
sol_0 = solutions(result_0)
HomotopyContinuation.write_parameters("Solution_0.txt", sol_0[1])
#sol_0 = HomotopyContinuation.read_parameters("Solution_0.txt")
=#

S_target_random1 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/S_random1.txt")
S_target_random2 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/S_random2.txt")
S_target_random3 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/S_random3.txt")

sol_0 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/Solution_0.txt")
sol_1a = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/Solution_1a.txt")
sol_2 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/Solution_2.txt")
sol_3 = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/Solution_3.txt")
sol_1b = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Monodromy/Solution_1b.txt")

L_tilde_0 =  L_start + sum(sol_0[i]*A_rand[i] for i in eachindex(A_rand)) 
L_tilde_1a = L_start + sum(sol_1a[i]*A_rand[i] for i in eachindex(A_rand)) 
L_tilde_1b = L_start + sum(sol_1b[i]*A_rand[i] for i in eachindex(A_rand)) 
L_tilde_2 = L_start + sum(sol_2[i]*A_rand[i] for i in eachindex(A_rand)) 
L_tilde_3 = L_start + sum(sol_3[i]*A_rand[i] for i in eachindex(A_rand)) 



verify_slicing(S_target_random3, L_tilde_3) 


minors(L_tilde_1a, 7)

a1 = det(L_tilde_1a[1:7,:])
b1 = det(L_tilde_1b[1:7,:])

a2 = det(L_tilde_1a[2:8,:])
b2 = det(L_tilde_1b[2:8,:])

a1/b1 == a2/b2 ## false => L_tilde_1a and L_tilde_1b are not the same up to GL(7)

