restart
notify = true
loadPackage "RealRoots"
load "AlgebraUtils.m2"

sliceComplementary = I -> (
    L = for i from 1 to dim I - 1 list random(1, ring I);
    J = I + ideal(L);
    return dehomogenize(J)
)

(k,n) = (2,6)
I = Grassmannian(k-1, n-1, CoefficientRing => QQ, Variable => p);

result = for i from 1 to 100 list (
    m = sliceComplementary I;
    f = quotMinPol(m,14);
    SturmCount(f)
)
tally result
