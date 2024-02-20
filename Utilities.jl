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
    R, t = PolynomialRing(K, varstring)
    vrs = t
    MS = MatrixSpace(R,k,m)
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

function poly_to_fp(F)
    n = Oscar.nvars(parent(F[1]))
    @var x_HC[1:n];
    F_HC = System([oscar_to_HC_Q(f, x_HC) for f in F])
    function F_fp(x)
        return HomotopyContinuation.ModelKit.evaluate(F_HC, x)
    end
    return F_fp
end


function Cayley_orth_to_skew(Q)
    #I = identity_matrix(QQ,nrows(Q))
    I = diagm([1 for i=1:7])
    Sk = inv(Q+I)*(Q-I)
    return Sk    
end


println("Include Utilities.jl done.")