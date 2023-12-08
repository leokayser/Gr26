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
Gamma[:,8:14]*transpose(Gamma[:,8:14])

R, var_x =  graded_polynomial_ring(QQ, ["x$i" for i=0:6]) 
#R = Q[x_i]
        
J = intersect([ideal(minors(hcat(matrix(var_x),Gamma[:,i]),2))  for i=1:14]);
I = intersect([ (ideal(minors(hcat(matrix(var_x),Gamma[:,i]),2))^2)  for i=1:14]);


betti(free_resolution(I))

####### problems
#=
gens3 = gens(truncate(J,3));
indx = findfirst(f -> (Oscar.total_degree(f)==3 ) , gens3) ; ## error 
gens3[indx]


is_degree_3(poly) = Oscar.total_degree(poly) ==3
index = findfirst(is_degree_3, gens3)
if index !== nothing
    println("First polynomial of degree 3: ", polynomials[index])
else
    println("No polynomial of degree 3 found.")
end
=#

B, _ = quo(R, I);
hilbert_function(B,3) #HF(R/I, 3)= 83
dim(collect(homogeneous_component(R, 3))[1]) # =84
# => J has 1 element of degree 3

minimal_betti_table(J)
##############


S, var_ax = graded_polynomial_ring(QQ, vcat(["a$i$j" for i=1:7 for j=1:15],["x$i" for i=0:6])) 
#S = Q[a_ij,x_i]

Sk = [ strictly_upper_triangular_matrix(var_ax[1+i*15:15+i*15]) - transpose(strictly_upper_triangular_matrix(var_ax[1+i*15:15+i*15])) for i=0:6  ]

A = sum( [(diagonal_matrix(var_ax[105+i],6))* (Sk[i]) for i=1:7])

#pf4 = pfaffians(A,4);
pf = pfaffian(A);
parent(pf)


#ι = hom(R,S, var_ax[106:112]) #inclusion from R to S
#G = pf-ι(g);
#base_ring(G) #QQ

Ta, var_a = graded_polynomial_ring(QQ, ["a$i$j" for i=1:7 for j=1:15] )  
##QQ[a_i_j]
Tax, var_ax = graded_polynomial_ring(Ta, ["x$i" for i=0:6])
##  QQ[a_i_j][x_i]

ι = hom(S,Tax,vcat(var_a, var_ax) ) 
Pf = ι(pf);
parent(Pf)

equations = []
for k=1:7
push!(equations, [ Oscar.evaluate(derivative(Pf,var_ax[k]), [Gamma[i,j] for i=1:nrows(Gamma)] ) for j=1:ncols(Gamma) ]) 
end

#Eq = hcat(equations...)
#length(equations[1])




#coeff = collect(Oscar.coefficients(ipf)); #  these are poly in a_ij
#parent(coeff[1])

## equation to set: coeff = 0 ( equations in Ta)


@var A[1:7,1:15] 
HC_vars = reshape(A[1:7,1:15] ,1,105)
coeff_HC =[ oscar_to_HC_Q(coeff[i],HC_vars) for i=1:length(coeff) ];

Eq = System(coeff_HC);
result = HomotopyContinuation.solve(Eq) #dim=21 underdetermined, makes sense





