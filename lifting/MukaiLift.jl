using Oscar
using LinearAlgebra
using HomotopyContinuation
#using AbstractAlgebra



function skew_matrix(a)
    U = strictly_upper_triangular_matrix(a)
    return U - transpose(U)
end

function randomSkmatrix(n)
    #S = matrix_space(QQ,n,n)
    v = [ rand(-8:8)//1 for i=1:n*(n-1)/2]
    return skew_matrix(v)
end

function randomSOmatrix(n)
    Sk = randomSkmatrix(n)
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


function max_independent_points(Γ)
    Γnew = Γ[:,1]
    ind = [1]
    i=2
    while (rank(Γnew)<=6)
        if rank(Γ[:,vcat(ind,[i])])>rank(Γnew)
            ind = push!(ind,i)
        end
        i = i+1
        Γnew = Γ[:,ind]
    end 
    return ind
end

function normalize_SApoints(Γ)
    ind = max_independent_points(Γ)
    Γ = hcat(Γ[:,ind],Γ[:,deleteat!([i for i=1:14],ind)])
    return inv(Γ[:,1:7])*Γ
end

function oscar_to_HC_Q(f,varstring)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(varstring.^exps[i]) for i = 1:length(cffs)])
end

function Cayley_skew_to_orth(Sk)
    I = identity_matrix(QQ,nrows(Sk))
    Q =  (I+Sk)*inv(I-Sk)
    return Q
end



function Cayley_skew_to_orth_coord(n)
    m = n*(n-1)÷ 2
    R, vrs = RationalFunctionField(QQ,["a_$i" for i=1:m])
    Sk = skew_matrix(vrs)
    I = identity_matrix(QQ,n)
    A = (I+Sk)
    B = I-Sk
    binv = [ (-1)^(i+j)*det(B[deleteat!(collect(1:n),i),deleteat!(collect(1:n),j)]) for i=1:n for j=1:n ]
    Binv = (1//(det(B))).*reshape( binv,n,n )
    Q = R.(A)*(Binv)
    #Jac = [ derivative(Q[i,j], vrs[k]) for i=1:n for j=1:n for k=1:length(vrs) ]
    return R,Q
end


 function Cayley_orth_to_skew(Q)
    I = identity_matrix(QQ,nrows(Q))
    Sk = inv(Q+I)*(Q-I)
    return Sk    
end

R,C = Cayley_skew_to_orth_coord(6)
R.(transpose(C))*C
C[1,1]

Γ = randomSApoints(14)

M = MatrixSpace(QQ,7,7)
A = M(reshape([rand(-10:10) for i=1:49],7,7))
G = hcat(A, A*(Γ[:,8:14]))
normalize_SApoints(G)== Γ # true

Sk = randomSkmatrix(6)
Q = Cayley_skew_to_orth(Sk)
Cayley_orth_to_skew(Q) == Sk #true

Q = randomSOmatrix(6)
Sk = Cayley_orth_to_skew(Q)
Cayley_skew_to_orth(Sk) == Q #true
