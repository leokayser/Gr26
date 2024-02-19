using Oscar
using LinearAlgebra
using HomotopyContinuation

println("Loading packages done!")

function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponents(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

function plückercoordinates(k,m,K)
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
    O = inv(Γ[:,1:7])*Γ[:,8:14]
    Γnorm = diagm([1 for i=1:7])
    for i = 1:7
        v = O[:,i]
        dot_prod = transpose(v)*v
        Γnorm = hcat(Γnorm, (1/sqrt(dot_prod)).*v)
        println(Γnorm)
    end
    return Γnorm
end


function Cayley_orth_to_skew(Q)
    #I = identity_matrix(QQ,nrows(Q))
    I = diagm([1 for i=1:7])
    Sk = inv(Q+I)*(Q-I)
    return Sk    
end


R,ϕ,vrs,M = plückercoordinates(2,6,QQ);

gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

A = reshape(gr_start_param[2:length(gr_start_param)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
B = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6
#B_normal_form = B*inv(B[1:7,:])


function poly_to_fp(F)
    n = Oscar.nvars(parent(F[1]))
    @var x_HC[1:n];
    F_HC = System([oscar_to_HC_Q(f, x_HC) for f in F])
    function F_fp(x)
        return HomotopyContinuation.ModelKit.evaluate(F_HC, x)
    end
    return F_fp
end


numerical_plücker = poly_to_fp(ϕ) 
Z = numerical_plücker.(gr_start_system) # 14 points in P^14

Γ = reshape(hcat([B\z for z in Z]...),7,14)


####################
Γ_norm = normalize_SApoints(Γ)
Γ_norm[:,8:14]*transpose(Γ_norm[:,8:14])


G = randomSApoints(14)
typeof(G)

M1 = reshape(randn(49),7,7)
D = diagm(randn(14))

M1*Float.(G)


Cayley_orth_to_skew(G[:,8:14])
S1 = Cayley_orth_to_skew(Γ_norm[:,8:14])

S1 + transpose(S1)
####################