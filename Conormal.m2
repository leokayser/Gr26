restart

n = 6
L = sort(subsets(1..n,2)) / (I -> 10*I_0 + I_1)
R = QQ[L/(I->p_I), L/(I->q_I), Degrees => toList( (#L:{0,1}) | (#L:{1,0}))]

P = genericSkewMatrix(R,p_12,n)
Q = genericSkewMatrix(R,q_12,n)

I = ideal(P*Q) + ideal(Q*P) + pfaffians(4,P) + pfaffians(6,Q);
betti mingens I
codim I, degree I
flatten entries (coefficients multidegree I)_1

M = monomialIdeal(I);
M == radical(M)
toString decompose M
#(decompose M)