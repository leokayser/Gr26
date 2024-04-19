using Oscar
using LinearAlgebra
using HomotopyContinuation

function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

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

function leadexp(f,w)
    exps = collect(Oscar.exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    exps[lm]
end

function toric_degen_poly_HC(f,w,vars_HC)
    terms = collect(Oscar.terms(f))
    fu = sum([ u^( transpose(w)*(leadexp(terms[i],w)-leadexp(f,w) )) * terms[i]  for i in eachindex(terms) ])
    fu_HC= oscar_to_HC_Q(fu, vars_HC) 
    return fu_HC
end


# input: a list of Oscar polynomials F 
# the polynomials in F are converted into HomotopyContinuation F_HC
# output: function evaluating the HomotopyContinuation polynomials in F_HC into a given input
function poly_to_fp(F)
    n = Oscar.nvars(parent(F[1]))
    @var x_HC[1:n];
    F_HC = System([oscar_to_HC_Q(f, x_HC) for f in F])
    function F_fp(x)
        return HomotopyContinuation.ModelKit.evaluate(F_HC, x)
    end
    return F_fp
end

#=
function Cayley_orth_to_skew(Q)
    #I = identity_matrix(QQ,nrows(Q))
    I = diagm([1 for i=1:7])
    Sk = inv(Q+I)*(Q-I)
    return Sk    
end

function Cayley_skew_to_orth(Sk)
    #I = identity_matrix(QQ,nrows(Sk))
    I = diagm([1 for i=1:7])
    Q =  (I+Sk)*inv(I-Sk)
    return Q
end
=#

function cayley_num(A)
    return (I+A)\(I-A) # = (I+A)^-1*(I-A)
end

function fake_Cayley_polynomial(n,vrs)
    #m = n*(n-1)÷ 2
    #R, vrs = RationalFunctionField(QQ,["a_$i" for i=1:m])
    Sk = skew_matrix(vrs)
    I = identity_matrix(QQ,n)
    numerator   = I-Sk
    denominator = I+Sk
    adjunct = reshape( [ (-1)^(i+j)*det(denominator[deleteat!(collect(1:n),i),deleteat!(collect(1:n),j)]) for i=1:n for j=1:n ], n,n)
    # no need to normalize once we are in P^6
    Q = Matrix(numerator)*adjunct
    return Q
end


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



function QQMatrix_to_ComplexF64(A)
    return ComplexF64.(Rational.(A))
end


function self_duality_control(Γ)
   B = [ Γ[i,1]*Γ[j,1] for i=1:7 for j=i:7]
   for h = 2:14
    C = [ Γ[i,h]*Γ[j,h] for i=1:7 for j=i:7]
    B = hcat(B,C)
   end
   λ = diagm(nullspace(B)[:,1])
   return λ
end

function normalize_SApoints(Γ)
    ind = max_independent_points(Γ)
    Γ = hcat(Γ[:,ind],Γ[:,deleteat!([i for i=1:14],ind)])
    λ = self_duality_control(Γ)
    λ = -im* λ
    λ_normalizing = diagm(vcat([sqrt(λ[i,i]) for i=1:7],[sqrt(-λ[i,i]) for i=8:14]))
    Γ_scaled = Γ * λ_normalizing
    Or = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    if rank(I+Or)<7
        λ = -im* λ
        λ_normalizing = diagm(vcat([sqrt(λ[i,i]) for i=1:7],[sqrt(-λ[i,i]) for i=8:14]))
        Γ_scaled = Γ * λ_normalizing
        Or = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    end
    Γnorm = hcat( diagm([1 for i=1:7]) , Or )
    return (Γnorm, Γ_scaled[:,1:7], λ_normalizing)
end

function proj_pts_eq(p,q)
    return (rank(hcat(p,q)) <= 1)
end

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

println("Include Utilities.jl done.")