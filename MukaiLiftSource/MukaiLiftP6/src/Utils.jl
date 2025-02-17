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

# Compute the Pluecker coordinates of Grassmannian Gr(2,m)
# -------------  Input:
# m         an integer
# K         a field
# -------------  Output:
# R         a polynomial ring
# ϕ         the pluecker coordinates
# vrs       the variables of the pluecker coordinates
# M         a matrix
function plückercoordinates(m,K)
    k = 2
    n = k*(m-k)
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
function poly_to_func(F)
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
function cayley(A)
    return (I+A)\(I-A) # = (I+A)^-1*(I-A)
end


# Compute skew normal form of A 
# -------------  Input:
# A            a matrix
# -------------  Output:
#              a matrix
function skew_to_SNF(A)
    return [I+A I-A]
end

# Compute a skew-symmetric matrix with entries in a
# -------------  Input:
# a            a vector
# -------------  Output:
#              a skew-symmetric matrix 
function vector_to_skew(a)
    U = strictly_upper_triangular_matrix(a)
    return U - transpose(U)
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



# Compute a 14x14 scaling matrix λ for a self-dual configuration Γ in P^6 verifying Γ*λ*Γ^T=0
# -------------  Input:
# Γ           a matrix whose columns represent self-dual points in P^6
# -------------  Output:
# λ           a matrix
function certify_selfdual(Γ)
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
# Γ_scaled[1:7]    the first 7 columns of a scaled form of Γ such that certify_selfdual(Γ_scaled)=I⊕(-I) 
# λ_normalizing    a scaling 14x14 matrix such that Γ_scaled= Γ*λ_normalizing
function orthogonal_normal_form(Γ)
    #ind = max_independent_points(Γ)
    #Γ = hcat(Γ[:,ind],Γ[:,deleteat!([i for i=1:14],ind)])
    λ_init = certify_selfdual(Γ)
    λ_normalizing = diagm(vcat([sqrt(Complex(λ_init[i,i])) for i=1:7],[sqrt(Complex(-λ_init[i,i])) for i=8:14]))
    Γ_scaled = Γ * λ_normalizing
    Or = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    #inv(I+Or)
    Γ_ONF = hcat( Diagonal([1 for i=1:7]) , -Or ) #!!!! "-" helps?
    A = Γ_scaled[:,1:7]
    λ = λ_normalizing * diagm(vcat(ones(7), -ones(7)))
    return (Γ_ONF, A, λ)
end

