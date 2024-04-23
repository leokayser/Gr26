include(pwd()*"/Utilities.jl")

function make_start()
    _,ϕ,_,_ = plückercoordinates(2,6,QQ);

    a1_start = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
    start_sol = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

    A = reshape(a1_start[2:length(a1_start)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
    L = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6


    numerical_plücker = poly_to_fp(ϕ) # Function that evaluates the polynomials in ϕ into the input
    Z = reshape(vcat(numerical_plücker.(start_sol)...),15,14)
    Γ = L\Z

    Γ_norm, A, λ = normalize_SApoints(Γ)

    norm(A*Γ_norm - Γ*λ, Inf)

    L_start = L*A

    norm(L_start*Γ_norm - Z*λ, Inf)

    Or = Γ_norm[:,8:14]

    S_start = cayley_num(Or)

    function simon(A)
        return [I+A I-A]
    end

    Gamma = simon(S_start)
    A2 = Gamma[:,1:7]
    norm(inv(A2)*Gamma - Γ_norm, Inf)

    L1 = L_start*inv(A2)
    norm(L1*Gamma - Z*λ, Inf)

    norm(cayley_num(S_start) - Or, Inf)
    skew_to_vector(S_start)

    return skew_to_vector(S_start), L1
end


