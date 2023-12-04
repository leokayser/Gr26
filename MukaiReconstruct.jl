using Oscar
using LinearAlgebra

function randomSOmatrix(n)
    #S = matrix_space(QQ,n,n)
    v = [ rand(-8:8)//1 for i=1:n*(n-1)/2]
    A = strictly_upper_triangular_matrix(v)
    Sk = A-transpose(A)
    I = identity_matrix(QQ,n)
    O = (I+Sk)^-1*(I-Sk)
    return O
end

function randomSApoints(s) # n is even, points are in P^(n/2)-1
    n = s÷2
    O = randomSOmatrix(n)
    pts = hcat(identity_matrix(QQ,n), O)
    return pts
end

Gamma = randomSApoints(14)

R, var =  graded_polynomial_ring(QQ, ["x$i" for i=0:6])

I = intersect([ideal(minors(hcat(matrix(var),Gamma[:,i]),2))  for i=1:14])

minimal_betti_table(I)