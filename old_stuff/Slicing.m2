restart
notify = true
loadPackage "RealRoots"
load "AlgebraUtils.m2"

sliceComplementary = I -> (
    L = for i from 1 to dim I - 1 list random(1, ring I);
    J = I + ideal(L);
    return dehomogenize(J)
)

(k,n) = (2,4)
I = Grassmannian(k-1, n-1, CoefficientRing => QQ, Variable => p);

result = for i from 1 to 100 list (
    m = sliceComplementary I;
    f = quotMinPol(m);
    SturmCount(f)
)
tally result

A = (1/2) * matrix{{0,1,0,0,0,0},{1,0,0,0,0,0},{0,0,0,1,0,0},{0,0,1,0,0,0},{0,0,0,0,0,1},{0,0,0,0,1,0}}
R = QQ[x_0..x_5]
X = matrix{{x_0..x_5}}
X*A*transpose(X)
S = QQ[p_0..p_14]
P = matrix{{p_0..p_14}}
(flatten entries(P*B*transpose(P)))_0