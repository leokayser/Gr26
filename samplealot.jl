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
    M = R.(hcat(M, reshape(t,k,m-k)))
    M = MS(M)
    ϕ = minors(M,k);
    return R,ϕ,vrs,M
end

function leadexp(f,w)
    exps = collect(Oscar.exponents(f))
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

Fu_HC = a*ϕu_HC #a are the linear forms cutting out a P^6 in P^14

C = System(Fu_HC, variables = vars_HC[1:8], parameters = [vars_HC[9];a[:]])

println("Defining system done.")

###########

#=
targ_par = [0;randn(ComplexF64, length(a[:]))]
torus_result = HomotopyContinuation.solve(C, target_parameters = targ_par )
torus_sol = solutions(torus_result);

targ_par_new = [1;randn(ComplexF64, length(a[:]))]
grass_result = HomotopyContinuation.solve(C, torus_sol; start_parameters = targ_par, target_parameters = targ_par_new)
grass_sol = solutions(grass_result)
HomotopyContinuation.write_parameters("Gr26_start_parameters.txt", targ_par_new)
HomotopyContinuation.write_solutions("Gr26_start_system.txt", grass_sol)
=#

############

   
gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

A = reshape(gr_start_param[2:length(gr_start_param)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
B = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6

Z = ϕ.gr_start_system

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
numerical_plücker.(gr_start_system) # 14 points in P^14


#=

println("Loading start parameters done!")

function uniform_Gr_point(a)
    return randn(length(a[:]))
end

N = parse(Int64, ARGS[1])
#N = 1000000
counter = Dict(2*i => 0 for i=0:7)

for i=1:N
    targ_par = [1; uniform_Gr_point(a)]
    grass_result = HomotopyContinuation.solve(C, gr_start_system; start_parameters = gr_start_param, target_parameters = targ_par)
    real_sols = length(HomotopyContinuation.real_solutions(grass_result))
    counter[real_sols] += 1
end

file = open(ARGS[2], "w")
for i in 0:7
    write(file, string(counter[2*i]),"\n")
end
close(file)

