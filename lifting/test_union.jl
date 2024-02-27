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