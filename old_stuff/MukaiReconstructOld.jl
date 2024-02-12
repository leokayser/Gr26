
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


function oscar_to_HC_Q(f,varstring)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(varstring.^exps[i]) for i = 1:length(cffs)])
end
 

Gamma = randomSApoints(14)



S, var_ax = polynomial_ring(QQ, vcat(["a_$i" for i=1:21],["x$i" for i=0:6])) 
#S = Q[a_ij,x_i]


function our_matrix_shape(i,vars)
    δ(x, y) = x==y ? 1 : 0
    lins = [transpose(rand(-10:10,3))*vars for i=1:5 ]
    v = S.([lins[1], δ(i,0), δ(i,1), vars[1], vars[2],
                    lins[2], δ(i,2), δ(i,3),  vars[3],
                            lins[3], δ(i,4),  δ(i,5),
                                    lins[4],  δ(i,6),
                                              lins[5]])
    return skew_matrix(v)
end


Sk = [ our_matrix_shape(i, var_ax[1+3*i:3+3*i]) for i=0:6  ]
#A0 = block_diagonal_matrix(QQ, [[0 1; -1 0],[0 1; -1 0],[0 0; 0 0]])
A =  sum([(diagonal_matrix(var_ax[21+i],6))* (Sk[i]) for i=1:7]);

#pf4 = pfaffians(A,4);
pf = pfaffian(A);
parent(pf)



Ta, var_a = polynomial_ring(QQ, ["a_$i" for i=1:21] )  
##QQ[a_i_j]
Tax, var_x = polynomial_ring(Ta, ["x$i" for i=0:6])
##  QQ[a_i_j][x_i]

ι = hom(S,Tax,vcat(var_a, var_x) ) 
Pf = ι(pf);
parent(Pf)
#collect(Oscar.terms(Pf))


equations = []
for k=1:7
push!(equations, [ Oscar.evaluate(derivative(Pf,var_x[k]), [Gamma[i,j] for i=1:nrows(Gamma)] ) for j=1:ncols(Gamma) ]) 
end


Eq = vcat(equations...);
length(Eq)

#I = ideal(Eq);

#a,b = Oscar.real_solutions(I);

### Convert to HomotopyContinuation
@var A[1:21] 
HC_vars = reshape(A[1:21] ,1,21)
Eq_HC =[ oscar_to_HC_Q(Eq[i],HC_vars) for i=1:length(Eq) ];
#for i=eachindex(Eq)
#   print(i);
#   oscar_to_HC_Q(Eq[i],HC_vars);
#end

result = HomotopyContinuation.solve( System(Eq_HC) ) 



