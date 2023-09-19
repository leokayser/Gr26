using Oscar
using LinearAlgebra
using HomotopyContinuation

println("Loading packages done!")

function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponent_vectors(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

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

function leadexp(f,w)
    exps = collect(exponents(f))
    weights = [dot(w,e) for e in exps]
    lm = argmin(weights)
    exps[lm]
end

k = 2
m = 6
R,ϕ,vrs,M = plückercoordinates(k,m,QQ);
T , variable = PolynomialRing(QQ, vcat(["t$i" for i=1:8], "u" ) )
u = variable[9]
ι = hom(R,T,gens(T)[1:8])
ϕ = ι.(ϕ)
w = [-3;0;-2;-1;-1;-2;0;-3;0];

terms = [ collect(Oscar.terms(ϕ[j])) for j in eachindex(ϕ) ]
ϕu = [sum([ u^( transpose(w)*(leadexp(terms[j][i],w)-leadexp(ϕ[j],w) )) * terms[j][i]  for i=1:length(terms[j]) ]) for j=1:length(ϕ)]
#l = 8
#Fu = [rand(-100:100,length(ϕu))'*ϕu for i = 1:l]

@var x[1:9] a[1:8,1:15]
vars_HC =  x[1:9] 
ϕu_HC= [oscar_to_HC_Q(ϕu[i], vars_HC) for i=1:length(ϕu)]

Fu_HC = a*ϕu_HC

C = System(Fu_HC, variables = vars_HC[1:8], parameters = [vars_HC[9];a[:]])

println("Defining system done.")


gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

println("Loading start parameters done!")

targ_par_int = [1; rand(-100:100, length(a[:]))]
@time grass_result = HomotopyContinuation.solve(C, gr_start_system; start_parameters = gr_start_param, target_parameters = targ_par_int)
HomotopyContinuation.real_solutions(grass_result)