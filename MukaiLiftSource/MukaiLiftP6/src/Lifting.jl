# Create a self-dual configuration with its embedding 
# -------------  Output:
# skew_to_vector(S_start)    the entries of the skew-symmetric matrix defining a configuration of self-dual points
# L1                         the embedding of the configuration in P14
function make_start()
    _,ϕ,_,_ = plückercoordinates(6,QQ);

    a1_start = HomotopyContinuation.read_parameters("MukaiLiftP6/src/Gr26_start_parameters.txt")
    start_sol = HomotopyContinuation.read_solutions("MukaiLiftP6/src/Gr26_start_system.txt")

    # A is the linear system whose solutions form a linear subspace L∼P^6 in P^14 
    A = reshape(a1_start[2:length(a1_start)], 8, 15) 
    L = LinearAlgebra.nullspace(A) 

    numerical_plücker = poly_to_func(ϕ) # Function that evaluates the polynomials in ϕ into the input
    Z = reshape(vcat(numerical_plücker.(start_sol)...),15,14)
    
    # Γ is the set of self-dual points embedded in L 
    Γ = L\Z
    Γ_ONF, A, λ = orthogonal_normal_form(Γ)
    
    @assert(norm(A*Γ_ONF - Γ*λ, Inf) < 1e-10)

    L_start = L*A
    
    #norm(L_start*Γ_ONF - Z*λ, Inf)

    Or = Γ_ONF[:,8:14]

    S_start = cayley(Or)
  
    Gamma = skew_to_SNF(S_start)
    A2 = Gamma[:,1:7]
    
    #norm(inv(A2)*Gamma - Γ_ONF, Inf)

    L1 = L_start*inv(A2) 

    #norm(L1*Gamma - Z*λ, Inf)
    #norm(cayley(S_start) - Or, Inf)
    
    skew_to_vector(S_start)

    return skew_to_vector(S_start), L1
end


# Compute the polynomial system and the start solution
# -------------  Output:
# parametrized_system     the system of equations that solve the Mukai lifting problem
# l_start                 a vector of zeros, the start solution 
function make_poly_system(S_start, L_start)

    @var l[1:69]
    @var s[1:21]
    @var x[1:7]
    @var q[1:15]

    S, s_oscar = polynomial_ring(QQ, ["s$i" for i=1:21])
    OscarSkewMat = vector_to_skew(s_oscar)
    Smat = [oscar_to_HC_Q(m, s) for m in OscarSkewMat]
    Γ = [I+Smat I-Smat]

    ######## To create different random matrices  ########
    #A_rand = [randn(ComplexF64,15,7) for _ in eachindex(l) ];
    #JLD.save("MukaiLiftP6/src/Random_matrices.jld", "data", A_rand )
    ######################################################

    A_rand = JLD.load("MukaiLiftP6/src/Random_matrices.jld")["data"];  

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

    return parametrized_system, l_start, A_rand
end