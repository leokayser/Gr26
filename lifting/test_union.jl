using HomotopyContinuation
include("union_system.jl")

@var x y z
F1 = System([x^2-1, y^2-2], variables=[x,y,z])
F2 = System([z^3-2], variables=[x,y,z])
F = [F1,F2]
sys = union_system(F)
size(sys)
a = [1,2,3]
sys(a)
res = zeros(ComplexF64, size(sys,1))
resJ = zeros(ComplexF64, size(sys))
HomotopyContinuation.ModelKit.evaluate_and_jacobian!(res,resJ,sys,a)
result = HomotopyContinuation.solve(sys)

@var x[1:2] a[1:2]
F1 = System([x[1]^2-a[1],x[1]-x[1]], variables = x, parameters = a)
F2 = System([x[1]*x[2]-a[1]+a[2]], variables = x, parameters = a)
F = union_system([F1,F2])
start_solutions = [[1, 1]]
p₁ = [1, 0]
p₀ = [2, 4]
HomotopyContinuation.solve(F, start_solutions; start_parameters=p₁, target_parameters=p₀)