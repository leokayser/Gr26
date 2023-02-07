using Oscar

using HomotopyContinuation

@var x[0:5]

#R = PolynomialRing(QQ, x)

#var(R)

G = grassmann_pluecker_ideal(2, 4)
L = [rand(-100:100,6)'*x for i = 1:4]

f=vcat(gens(G),L)

typeof(gen)

variables(f)
F = System(f; variables = x)  