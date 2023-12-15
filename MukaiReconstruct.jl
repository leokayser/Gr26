
using Oscar
using LinearAlgebra
using HomotopyContinuation
#using AbstractAlgebra

function randomSkmatrix(n)
    #S = matrix_space(QQ,n,n)
    v = [ rand(-8:8)//1 for i=1:n*(n-1)/2]
    A = strictly_upper_triangular_matrix(v)
    Sk = A-transpose(A)
    return Sk
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


function oscar_to_HC_Q(f,varstring)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(varstring.^exps[i]) for i = 1:length(cffs)])
end

Gamma = randomSApoints(14)



S, var_ax = graded_polynomial_ring(QQ, vcat(["a_($i,$j)" for i=1:6 for j=1:15],["x$i" for i=0:6])) 
#S = Q[a_ij,x_i]

Sk = [ strictly_upper_triangular_matrix(var_ax[1+i*15:15+i*15]) - transpose(strictly_upper_triangular_matrix(var_ax[1+i*15:15+i*15])) for i=0:5  ]
A0 = block_diagonal_matrix(QQ, [[0 1; -1 0],[0 1; -1 0],[0 0; 0 0]])
A = diagonal_matrix(var_ax[91],6)*A0 + sum( [(diagonal_matrix(var_ax[91+i],6))* (Sk[i]) for i=1:6]);

#pf4 = pfaffians(A,4);
pf = pfaffian(A);
parent(pf)



Ta, var_a = graded_polynomial_ring(QQ, ["a$i;$j" for i=1:6 for j=1:15] )  
##QQ[a_i_j]
Tax, var_ax = graded_polynomial_ring(Ta, ["x$i" for i=0:6])
##  QQ[a_i_j][x_i]

ι = hom(S,Tax,vcat(var_a, var_ax) ) 
Pf = ι(pf);
parent(Pf)
#collect(Oscar.terms(Pf))


equations = []
for k=1:7
push!(equations, [ Oscar.evaluate(derivative(Pf,var_ax[k]), [Gamma[i,j] for i=1:nrows(Gamma)] ) for j=1:ncols(Gamma) ]) 
end

Eq = hcat(equations...);
length(Eq)


### Convert to HomotopyContinuation
@var A[1:6,1:15] 
HC_vars = reshape(A[1:6,1:15] ,1,90)
Eq_HC =[ oscar_to_HC_Q(Eq[i],HC_vars) for i=1:length(Eq) ];
for i=eachindex(Eq)
   print(i);
   oscar_to_HC_Q(Eq[i],HC_vars);
end

oscar_to_HC_Q(Eq[2],HC_vars)

result = HomotopyContinuation.solve( System(Eq_HC) ) #dim=21 underdetermined, makes sense





