module KhovanskiiSolving

# Authors: Barbara Betti, Marta Panizzut and Simon Telen
# Data: May 22, 2023
# Short description: This code accompanies the paper
# [BPT] Barbara Betti, Marta Panizzut and Simon Telen, "Solving equations on parameterized varieties", arXiv:


using Oscar
using LinearAlgebra

export leadterm,
leadexp,
leadmon,
get_basis,
getHF,
expand_Mac,
get_KM,
max_independent_rows,
get_multiplication_matrices,
get_solutions,
get_real_solutions,
get_commuting_matrices,
solve_Khovanskii,
reduceMac,
HF_Gr,
deg_Gr,
plückercoordinates,
dimSchubertVariety,
codimSchubertVariety,
equationsSchubertVariety

# Compute the leading term of f with respect to the weight vector w.
function leadterm(f,w)
    exps = collect(exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    collect(terms(f))[lm]
end

# Compute the leading exponent of f with respect to the weight vector w.
function leadexp(f,w)
    exps = collect(exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    exps[lm]
end

# Compute the leading monomial of f with respect to the weight vector w.
function leadmon(f,w)
    exps = collect(exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    collect(Oscar.monomials(f))[lm]
end

#From now ϕ is a Khovanskii basis in the variables in vars. The leading exponents are in leadexps.

# Compute a basis for K[ϕ]_d, the degree d part of the algebra K[ϕ].
# The first output is the basis, the second the vector with the leading exponents of the elements in the basis
function get_basis(ϕ,d,K,vars,leadexps)
    vars_t = ["t$i" for i = 1:length(ϕ)]
    S, t = PolynomialRing(K, vars_t)
    allmons = collect(Oscar.monomials(sum(t)^d)) #all monomials in t of degree d
    allexps = collect(Oscar.exponents(sum(t)^d)) #all exponents of monomials in t of degree d
    latticepoint = [ allexps[i]'*leadexps    for i=1:length(allexps)  ] #all d-sum of elements in leadexps
    mons_in_dP = [prod(vars.^lp) for lp in latticepoint]  #monomials with exponent in latticepoint
    mons_in_dP_unique = unique(mons_in_dP)
    leadmons = [prod(vars.^e) for e in leadexps]
    aux = [findfirst(mon-> Oscar.evaluate(mon,leadmons) == mymon, allmons) for mymon in mons_in_dP_unique]
    [prod(ϕ.^e) for e in allexps[aux]], allexps[aux] #allexps[aux][i]:= allexps[aux[i]]
end

# Compute the Hilbert Function of K[ϕ] in degree d
function getHF(ϕ,d,K,vars,leadexps)
    return length(get_basis(ϕ,d,K,vars,leadexps)[1])
end

# Compute the Khovanskii-Macaulay matrix w.r.t the polynomials in F in degree d.
# The optional input are
function expand_Mac(F,shifts,exps,ϕ,K,vars; max_iter = 200, random_range = 50)
    L = length(exps)
    MSl = MatrixSpace(K,L,L); MSr = MatrixSpace(K,L,sum(length.(shifts)))
    Vdm = zero(MSl); rhs = zero(MSr)
    success = false
    iter = 1
    c = []
    while success == false && iter < max_iter
        if K==QQ
            pts = [rand(-random_range:random_range,length(vars)).//rand(1:random_range,length(vars) ) for i = 1:L] #vector with L elements
        else
            pts = [rand(0:characteristic(K)-1,length(vars)) for i = 1:L] #vector with L elements
        end
        for i = 1:L
            hval = [Oscar.evaluate(hh,pts[i]) for hh in ϕ] #vector with elements of polys evaluated in pts(i), length(hval)=length(polys)
            Vdm[i,:] = [prod(hval.^e) for e in exps] #prod([a..n].^[a'..n'])=a^a'\cdot..\cdot n^n'
            rhs[i,:] = [Oscar.evaluate(F[j],pts[i])*prod(hval.^shift) for j = 1:length(F) for shift in shifts[j]] # f(pts[i])* hval.^shift , this is a scalar
        end
        if rank(Vdm) == L
            success = true
        else
            iter +=1
        end
    end
    if success

        return transpose(solve(Vdm,rhs))
    else
        println("couldn't make an invertible matrix ")
    end
end


# Compute the Khovanskii-Macaulay matrix of the polynomials in F in degree d.
# degs_F is the vector with the degree of the polynomials in F in the Khovanskii-basis ϕ
# The optional input reduce is false by default. If reduce=true the Khovanskii-Macaulay matrix is reduced to a full rank submatrix
function get_KM(F,d,degs_F,ϕ,K,vars,leadexps; reduce=false )
    shifts= [get_basis(ϕ,d-degs_F[i],K,vars,leadexps)[2] for i = 1:length(F)]
    exps = get_basis(ϕ,d,K,vars,leadexps)[2]
    Mac = expand_Mac(F,shifts,exps,ϕ,K,vars)
    if reduce == true
        Mac=reduceMac(Mac,K)
    end
    return Mac
end


# The function selects a maximal set of independent rows of a matix M in ascending order.
# The first output is the corresponding full-rank submatix. The second output is the vector with the indices of the selected rows.
function max_independent_rows(M; K = QQ, cautious = false)
    start_ind = findfirst(j-> M[j,:] != [0 for i = 1:size(M,1)], 1:size(M,1)) #index of the first not zero row
    if !isempty(start_ind)
        inds = [start_ind]
        MM = M[start_ind,:]
        r = 1
        i = start_ind + 1
        while r < minimum(size(M)) && i <= size(M,1)
            MMnew = vcat(MM,M[i,:])
            if K == QQ && !cautious
                MMnew_int = maximum(denominator.(MMnew)).*numerator.(MMnew)
                MS = MatrixSpace(GF(508549),size(MMnew_int)...)
                MMnew_int = MS(MMnew_int)
            else
                MMnew_int = MMnew
            end
            rk = rank(MMnew_int)
            if rk > r
                r = rk
                MM = MMnew
                push!(inds,i)
            end
            i += 1
        end
        if rank(M) == rank(MM)
            return MM, inds
        else
            return max_independent_rows(M; K=K, cautious=true)
        end
    else
        return [], []
    end
end


# The first output is the vector with the multiplication matrices whose eigenvalues are the solutions to the system.
# The second output is the vector c that give the linear combination of the polynomials in ϕ used to compute the solutions
function get_multiplication_matrices(Mac,d,ϕ,K,vars,leadexps)
    N = nullspace(Mac)[2];
    Ns = []
    basisdminus1 = get_basis(ϕ,d-1,K,vars,leadexps)[1]
    exps = get_basis(ϕ,d,K,vars,leadexps)[2]
    for i = 1:length(ϕ)
        idmtx = Matrix(I,length(ϕ),length(ϕ))
        W = expand_Mac(basisdminus1,[[idmtx[i,:]] for j = 1:length(basisdminus1)],exps,ϕ,K,vars)
        #@show(size(W),size(N))
        push!(Ns,W*N)
    end
    c = rand(-100:100,length(ϕ));
    Nstar = sum([c[i]*Ns[i] for i = 1:length(ϕ)]);
    Mulstar, inds = max_independent_rows(Nstar; K = K);
    Mul = [solve(Mulstar,NN[inds,:]) for NN in Ns];
    return Mul, c
end


# The output are the parameterized solutions of the system when the field is Q
function get_solutions(Mul)
    Mul_rat = [convert.(Rational{BigInt},MM) for MM in Mul];
    Mul_float = [Float64.(MM) for MM in Mul_rat];
    Mul_rand = sum([randn()*Mf for Mf in Mul_float])
    eigenobj = eigen(Mul_rand)
    V = eigenobj.vectors
    sols = hcat([diag(inv(V)*Mf*V) for Mf in Mul_float]...)
    return sols
end

# The output are the real parameterized solutions of the system when the field is Q.
# The optional input is the value of the tolerance, that is 1e-10 by default.
# A solution is considered real when every coordinate has the ratio between imaginary part and norm that is less than the tolerance threshold
function get_real_solutions(Mul; tol = 1e-10)
    Mul_rat = [convert.(Rational{BigInt},MM) for MM in Mul];
    Mul_float = [Float64.(MM) for MM in Mul_rat];
    Mul_rand = sum([randn()*Mf for Mf in Mul_float])
    eigenobj = eigen(Mul_rand)
    V = eigenobj.vectors
    sols = hcat([diag(inv(V)*Mf*V) for Mf in Mul_float]...)
    real_inds = []
    for i = 1:size(sols,1)
        if norm(imag(sols[i,:]))/norm(sols[i,:]) < tol
            push!(real_inds,i)
        end
    end
    return real.(sols[real_inds,:])
end

# This function computes directly the multiplication matrices. The input are the list of polynomials F,
# the degree in which compute the Khovanskii-Macaulay matrix, the degree of polynomials in F, the basis ϕ, the field K, the variables and the leading exponents of ϕ
# The output is the vector with the multiplication matrices
function get_commuting_matrices(F,dreg,degs_F,ϕ,K,vars,leadexps)
    Mac = get_KM(F,dreg,degs_F,ϕ,K,vars,leadexps)
    Mul, c = get_multiplication_matrices(Mac,dreg,ϕ,K,vars,leadexps)
    return Mul
end

# This
function solve_Khovanskii(F,dreg,degs_F,ϕ,vars,leadexps)
    Mac = get_KM(F,dreg,degs_F,ϕ,QQ,vars,leadexps)
    Mac = reduceMac(Mac,QQ)
    Mul, c = get_multiplication_matrices(Mac,dreg,ϕ,QQ,vars,leadexps)
    return get_solutions(Mul)
end


function reduceMac(Mac,K)
    MSa=MatrixSpace(K,rank(Mac),size(Mac)[1]);
    T=MSa(rand(-100:100,rank(Mac),size(Mac)[1]));
    MacNew=T*Mac;
    return MacNew;
end


####### Schubert calculus

# Hilbert function of Grassmannian Gr(k,m)
function HF_Gr(k,m,t)
     prefac=prod([factorial(i) for i=1:k-1])//(prod([factorial(m-i) for i=1:k]))
     prefac*(prod([prod([t+i+l for l=0:m-k-1]) for i=1:k]))
end

# Degree of Grassmannian Gr(k,m)
function deg_Gr(k,m)
    D=k*(m-k)
    prefac=prod([factorial(i) for i=1:k-1])//(prod([factorial(m-i) for i=1:k]))
    deg=factorial(D)*prefac
end

# It returns the plucker relations of the Grassmannian Gr(k,m),
# i.e. kxk minors of matrix with identity on the left and variables on the right( I_k | t_i )
function plückercoordinates(k,m,K)
    n=k*(m-k)
    varstring = ["t$i" for i =1:n];
    R, t = PolynomialRing(K, varstring)
    vrs = t
    MS = MatrixSpace(R,k,m)
    M = Int64.(diagm(ones(k))).+0*t[1]
    M = R.(hcat(M, reshape(t,k,m-k)))
    M = MS(M)
    ϕ = minors(M,k);
    return R,ϕ,vrs,M
end

# Compute the dimension of the Schubert variety in Gr(k,m) defined by the Schubert condition α
function dimSchubertVariety(α,k)
    return (sum(α))-(k*(k+1))/2
end

# Compute the codimension of the Schubert variety in Gr(k,m) defined by the Schubert condition α
function codimSchubertVariety(α,k,m)
    n=k*(m-k)
    return n-dimSchubertVariety(α,k)
end

# Compute the defining equations of the Schubert variety defined by the Schubert conditions in A=[α_1,..,α_c].
# Returns the flags, the equations,
function equationsSchubertVariety(A,k,m,K)
    R,ϕ,vrs,M = plückercoordinates(k,m,K)
    MS1 = MatrixSpace(R,m,m);
    Flags = [ MS1( reshape(rand(-100:100,m^2),m,m) ) for i=1:length(A) ]
    F=[]
    for j=1:length(A)
        for i=1:length(A[j])
        if k+1+A[j][i]-i <=minimum([m,k+A[j][i]])
            MSi=MatrixSpace(R,k+A[j][i],m)
        push!(F, minors( MSi( vcat(M,Flags[j][[i for i=1:A[j][i]],:])), k+1+A[j][i]-i ) )
        end
        end
    end
    F=unique(vcat(F...));
    w = [-1 for i=1:length(vrs)];
    leadexps = [leadexp(hh,w) for hh in ϕ];
    d=1
    Mac = get_KM(F,d,[1 for i=1:length(F)],ϕ,K,vrs,leadexps);
    MSa=MatrixSpace(R,rank(Mac),size(Mac)[1])
    T=MSa(rand(-100:100,rank(Mac),size(Mac)[1]));
    F=T*F;
    degs_F=[1 for i=1:length(F)];
    return Flags,F,leadexps,degs_F
end


end # module KhovanskiiSolving
