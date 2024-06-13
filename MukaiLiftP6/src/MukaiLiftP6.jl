module MukaiLiftP6

# Authors: Barbara Betti and Leonie Kayser
# Data: Jun 13, 2024
# Short description: This code accompanies the paper
# [BK24] Barbara Betti and Leonie Kayser, "Mukai lifting of self-dual points in P6", arXiv: https://arxiv.org/abs/2406.02734

using Oscar
using LinearAlgebra
using HomotopyContinuation
using DelimitedFiles
using JLD


export leadexp,
plückercoordinates,
oscar_to_HC_Q,
toric_degen_poly_HC,
poly_to_fp,
cayley_num,
skew_normal_form,
skew_matrix,
randomSkmatrix,
randomSOmatrix,
randomSApoints,
max_independent_points,
QQMatrix_to_ComplexF64,
self_duality_control,
normalize_SApoints,
proj_pts_eq,
skew_to_vector,
make_start,
make_poly_system



# Compute the leading exponent of f with respect to the weight vector w.
# -------------  Input:
# f          a polynomial
# w          a weight vector
# -------------  Output:
# leading exponent of f 
function leadexp(f,w)
    exps = collect(Oscar.exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    exps[lm]
end

# Compute the Pluecker coordinates of Grassmannian Gr(k,m)
# -------------  Input:
# k         an integer
# m         an integer
# K         a field
# -------------  Output:
# R         a polynomial ring
# ϕ         the pluecker coordinates
# vrs       the variables of the pluecker coordinates
# M         a matrix
function plückercoordinates(k,m,K) # only works for k=2
    n=k*(m-k)
    varstring = ["t$i" for i =1:n];
    R, t = polynomial_ring(K, varstring)
    vrs = t
    MS = matrix_space(R,k,m)
    M = Int64.(diagm(ones(k))).+0*t[1]
    Q = transpose(reshape(t,m-k,k))
    M = R.(hcat(M, Q))
    M = MS(M)
    ϕ =[]
    for i=1:m-1
        for j=i+1:m
            push!(ϕ, det(M[:,[i,j]]))
        end
    end
    return R,ϕ,vrs,M
end

# Convert a Oscar polynomial in an expression HomotopyContinuation 
# -------------  Input:
# f          a polynomial ::QQPolyRingElem
# vars       HC variables ::Vector{Variable}
# -------------  Output:
#            HC polynoial ::expression 
function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

# Compute the HomotopyContinuation polynomial in an extra variables 
# u that degenerateds to the initial form of f for u->0
# -------------  Input:
# f          a polynomial ::QQPolyRingElem
# w          a weight vector
# vars_HC    HC variables ::Vector{Variable}
# u          parameter of the toric degeneration ::Variable
# -------------  Output:
# fu_HC      a HC polynomial ::expression 
function toric_degen_poly_HC(f,w,vars_HC,u)
    terms = collect(Oscar.terms(f))
    fu = sum([ u^( transpose(w)*(leadexp(terms[i],w)-leadexp(f,w) )) * terms[i]  for i in eachindex(terms) ])
    fu_HC= oscar_to_HC_Q(fu, vars_HC) 
    return fu_HC
end

# Compute the evaluation function of the polynomial in F in HomotopyContinuation
# -------------  Input:
# F          a vector of Oscar polynomials ::Vector{QQPolyRingElem}
# -------------  Output:
# F_fp        a function 
function poly_to_fp(F)
    n = Oscar.nvars(parent(F[1]))
    @var x_HC[1:n];
    F_HC = System([oscar_to_HC_Q(f, x_HC) for f in F])
    function F_fp(x)
        return HomotopyContinuation.ModelKit.evaluate(F_HC, x)
    end
    return F_fp
end

# Compute numerically the Cayley transform of A 
# -------------  Input:
# A            a skew-symmetric or orthogonal matrix
# -------------  Output:
#              an orthogonal or skew-symmetric matrix 
function cayley_num(A)
    return (I+A)\(I-A) # = (I+A)^-1*(I-A)
end


# Compute skew normal form of A 
# -------------  Input:
# A            a matrix
# -------------  Output:
#              a matrix
function skew_normal_form(A)
    return [I+A I-A]
end

# Compute a skew-symmetric matrix with entries in a
# -------------  Input:
# a            a vector
# -------------  Output:
#              a skew-symmetric matrix 
function skew_matrix(a)
    U = strictly_upper_triangular_matrix(a)
    return U - transpose(U)
end

# Compute a random skew-symmetric matrix of size nxn
# -------------  Input:
# n            an integer
# -------------  Output:
#              a nxn skew-symmetric matrix 
function randomSkmatrix(n)
    v = [ rand(-8:8)//1 for i=1:n*(n-1)/2]
    return skew_matrix(v)
end

# Compute a random special orthogonal matrix of size nxn
# -------------  Input:
# n            an integer
# -------------  Output:
#              a nxn special orthogonal matrix 
function randomSOmatrix(n)
    Sk = randomSkmatrix(n)
    I = identity_matrix(QQ,n)
    O = (I+Sk)^-1*(I-Sk)
    return O
end

# Compute a random configuration of s self-dual points
# -------------  Input:
# s            an even integer
# -------------  Output:
# pts          a s/2 x s matrix whose columns represent points in P^(n/2)-1
function randomSApoints(s) 
    n = s÷2
    O = randomSOmatrix(n)
    pts = hcat(identity_matrix(QQ,n), O)
    return pts
end


# Compute a maximal independent set of columns of a matrix Γ
# -------------  Input:
# Γ            a matrix
# -------------  Output:
# ind          a vector of integers indicating the indexes of the columns
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


# Convert the rational entries of a vector or a matrix to complex numbers
# -------------  Input:
# A            a vector or a matrix
# -------------  Output:
#              the matrix A with entries over C
function QQMatrix_to_ComplexF64(A)
    return ComplexF64.(Rational.(A))
end


# Compute a 14x14 scaling matrix λ for a self-dual configuration Γ in P^6 verifying Γ*λ*Γ^T=0
# -------------  Input:
# Γ           a matrix whose columns represent self-dual points in P^6
# -------------  Output:
# λ           a matrix
function self_duality_control(Γ)
   B = [ Γ[i,1]*Γ[j,1] for i=1:7 for j=i:7]
   for h = 2:14
    C = [ Γ[i,h]*Γ[j,h] for i=1:7 for j=i:7]
    B = hcat(B,C)
   end
   #print(nullspace(B))
   λ = diagm(nullspace(B)[:,1])
   return λ
end

# Compute the orthogonal normal form for a self-dual configuration Γ in P^6
# -------------  Input:
# Γ           a matrix whose columns represent self-dual points in P^6
# -------------  Output:
# Γ_norm           an orthogonal normal form of Γ   
# Γ_scaled[1:7]    the first 7 columns of a scaled form of Γ such that self_duality_control(Γ_scaled)=I⊕(-I) 
# λ_normalizing    a scaling 14x14 matrix such that Γ_scaled= Γ*λ_normalizing
function normalize_SApoints(Γ)
    ind = max_independent_points(Γ)
    Γ = hcat(Γ[:,ind],Γ[:,deleteat!([i for i=1:14],ind)])
    λ = self_duality_control(Γ)
    #λ = -im* λ # WHY IS THIS HERE???
    λ_normalizing = diagm(vcat([sqrt(Complex(λ[i,i])) for i=1:7],[sqrt(Complex(-λ[i,i])) for i=8:14]))
    Γ_scaled = Γ * λ_normalizing
    Or = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    #=if rank(I+Or)<7 # TODO: rework
        λ = -im* λ
        λ_normalizing = diagm(vcat([sqrt(λ[i,i]) for i=1:7],[sqrt(-λ[i,i]) for i=8:14]))
        Γ_scaled = Γ * λ_normalizing
        Or = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    end=#
    Γ_ONF = hcat( Diagonal([1 for i=1:7]) , Or )
    return (Γ_ONF, Γ_scaled[:,1:7], λ_normalizing)
end

# Control if two points are projectively equivalent
# -------------  Input:
# p           a vector
# q           a vector
# -------------  Output:
#           true if p and q are linearly dependent
function proj_pts_eq(p,q)
    return (rank(hcat(p,q)) <= 1)
end

# Compute the vector of entries of a skew-symmetric matrix
# -------------  Input:
# A            a nxn skew-symmetric matrix
# -------------  Output:
# v            a vector of length n
function skew_to_vector(A)
    n = ncols(A)
    v = []
    for i=1:n
        for j=(i+1):n
            push!(v, A[i,j])
        end
    end
    return v
end


# Create a self-dual configuration with its embedding 
# -------------  Output:
# skew_to_vector(S_start)    the entries of the skew-symmetric matrix defining a configuration of self-dual points
# L1                         the embedding of the configuration in P14
function make_start()
    _,ϕ,_,_ = plückercoordinates(2,6,QQ);

    a1_start = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
    start_sol = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

    # A is the linear system whose solutions form a linear subspace L∼P^6 in P^14 
    A = reshape(a1_start[2:length(a1_start)], 8, 15) 
    L = LinearAlgebra.nullspace(A) 

    numerical_plücker = poly_to_fp(ϕ) # Function that evaluates the polynomials in ϕ into the input
    Z = reshape(vcat(numerical_plücker.(start_sol)...),15,14)
    
    # Γ is the set of self-dual points embedded in L 
    Γ = L\Z
    Γ_ONF, A, λ = normalize_SApoints(Γ)
    #norm(A*Γ_ONF - Γ*λ, Inf) 

    L_start = L*A
    #norm(L_start*Γ_ONF - Z*λ, Inf)

    Or = Γ_ONF[:,8:14]

    S_start = cayley_num(Or)
  
    Gamma = skew_normal_form(S_start)
    A2 = Gamma[:,1:7]
    
    #norm(inv(A2)*Gamma - Γ_ONF, Inf)

    L1 = L_start*inv(A2) 

    #norm(L1*Gamma - Z*λ, Inf)
    #norm(cayley_num(S_start) - Or, Inf)
    
    skew_to_vector(S_start)

    return skew_to_vector(S_start), L1
end


# Compute the polynomial system and the start solution
# -------------  Output:
# parametrized_system     the system of equations that solve the Mukai lifting problem
# l_start                 a vector of zeros, the start solution 
function make_poly_system(S_start, L_start)

    @var l[1:69]
    @var s[1:21]
    @var x[1:7]
    @var q[1:15]

    S, s_oscar = polynomial_ring(QQ, ["s$i" for i=1:21])
    OscarSkewMat = skew_matrix(s_oscar)
    Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]
    Γ = [I+Smat I-Smat]

    A_rand = JLD.load("Random_matrices.jld")["data"];  

    L_l = L_start + sum(l[i]*A_rand[i] for i in eachindex(l));
    L_sys = System(L_l*x, variables=[x;l]);

    plück_oscar = gens( grassmann_pluecker_ideal(2,6))
    plück_sys = System([oscar_to_HC_Q(plück_oscar[i], q) for i=1:15], variables=q)

    Q = plück_sys(expressions(L_sys));
    Q_sys = System(Q, variables=[x;l]);
    equations = vcat([Q_sys([Γ[:,i];l]) for i in 1:14]...);

    parametrized_system = System(equations, variables=l, parameters=s);
    # norm( parametrized_system(l_start,S_start), Inf)   

    l_start = zeros(ComplexF64,length(l))

    return parametrized_system, l_start
end

end # module MukaiLiftP6
