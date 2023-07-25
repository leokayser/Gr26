using Oscar
using HomotopyContinuation
using StatsBase

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

function deg_Gr(k,m)
    D=k*(m-k)
    prefac=prod([factorial(i) for i=1:k-1])//(prod([factorial(m-i) for i=1:k]))
    deg=factorial(D)*prefac
end


k = 2
m = 6
f = random_linear_section(grassmann_pluecker_ideal(k,m), k*(m-k))
HC_solve_oscar(f)