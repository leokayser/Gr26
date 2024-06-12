include("Utilities.jl")

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

S_start, L_start = make_start();

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


######## Solving for a random configuration ########

# S_target = randn(ComplexF64,21)
# @time result = HomotopyContinuation.solve(R, l_start; start_parameters=S_start, target_parameters=S_target)

####################################################