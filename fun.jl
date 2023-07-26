using Oscar
using HomotopyContinuation
using KhovanskiiSolving
#using AbstractAlgebra

function oscar_to_HC_Q(f,vars)
    cffs = convert.(Rational{Int64},collect(Oscar.coefficients(f)))
    exps = collect(Oscar.exponent_vectors(f))
    sum([cffs[i]*prod(vars.^exps[i]) for i = 1:length(cffs)])
end

function random_linear_section(I,l)
    oscar_vars = gens(base_ring(I))
    L = [rand(-10:10,length(oscar_vars))'*oscar_vars for i = 1:l]
    return vcat(gens(I),L)
end

function HC_solve_oscar(f)
    N = length(gens(base_ring(ideal([f[1]]))))
    @var x[0:N-1]
    hc_f = [oscar_to_HC_Q(fi, x) for fi in f]
    F = System(hc_f; variables = x)  
    return HomotopyContinuation.solve(F)
end

function pluecker_embedding(mat)
    k = size(mat)[1]
    m = size(mat)[2]
    S = MatrixSpace(ZZ, k, m)
    return minors(S(mat), k)
end

function linear_forms(pts,ϕ)
    r = length(pts)
    N = length(pts[1])
    R = ϕ[1].parent
    mat = reshape(vcat(pts...), N,r)
    S = MatrixSpace(R,N,r+1)
    matx = S(hcat(mat, ϕ))
    lins = [det(matx[[1:r;row],:]) for row in (r+1):N]
    return lins
end

realpts = [pluecker_embedding(rand(-10:10,2,6)) for i = 1:7]

transpose(reshape( vcat(realpts...), 15,7))

k = 2
m = 6
f = random_linear_section(grassmann_pluecker_ideal(k,m), k*(m-k))
HC_solve_oscar(f)


###############################################
k = 2
m = 6
R,ϕ,vrs,M= plückercoordinates(k,m,QQ);
R
T , variable = PolynomialRing(QQ, vcat(["t$i" for i=1:8], "u" ) )
u = variable[9]
ι = hom(R,T,gens(T)[1:8])
ϕ = ι.(ϕ)
w = [-4;-1;-3;-2;-2;-3;-1;-4;0];

terms = [ collect(Oscar.terms(ϕ[j])) for j in eachindex(ϕ) ]
ϕu = [sum([ u^( transpose(w)*(leadexp(terms[j][i],w)-leadexp(ϕ[j],w) )) * terms[j][i]  for i=1:length(terms[j]) ]) for j=1:length(ϕ)]

l = 8
Fu = [rand(-10:10,length(ϕu))'*ϕu for i = 1:l];

 @var x[1:9]
vars_HC =  x[1:9]
Fu_HC = [oscar_to_HC_Q(F[j],vars_HC) for j in eachindex(Fu)];

C = System(Fu_HC, variables = vars_HC[1:8], parameters = vars_HC[9:9]);

ev = hom(T,T,vcat(gens(T)[1:8],[0]))
F0 = ev.(Fu)

torus_result = HomotopyContinuation.solve(C, target_parameters = [0])
torus_sol = solutions(torus_result);
HomotopyContinuation.solve(C, torus_sol; start_parameters = [0], target_parameters = [1])
