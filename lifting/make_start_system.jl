include(pwd()*"/Utilities.jl")

function make_start()
    _,ϕ,_,_ = plückercoordinates(2,6,QQ);

    gr_start_param = HomotopyContinuation.read_parameters("Gr26_start_parameters.txt")
    gr_start_system = HomotopyContinuation.read_solutions("Gr26_start_system.txt")

    A = reshape(gr_start_param[2:length(gr_start_param)], 8, 15) # This is the linear system whose solution are the Λ in P^14 
    φ = LinearAlgebra.nullspace(A) # each column is an element of the basis of the space  Λ in P^14. The basis determines an isomorphism with P^6
    #B_normal_form = B*inv(B[1:7,:])


    numerical_plücker = poly_to_fp(ϕ) # Function that evaluates the polynomials in ϕ into the input
    Z = numerical_plücker.(gr_start_system) # 14 points in P^14
    Z_mat = reshape(vcat(Z...),15,14)

    Γ = reshape(hcat([φ\z for z in Z]...),7,14)

    Γ_norm, A, λ = normalize_SApoints(Γ)

    norm(A*Γ_norm - Γ*λ, Inf)

    φ_start = φ*A

    norm(φ_start*Γ_norm - Z_mat*λ, Inf)

    #proj_pts_eq((φ_start*Γ_norm)[:,8],Z[8])

    O = Γ_norm[:,8:14]

    S_start = cayley_num(O)
    norm(cayley_num(S_start) - O, Inf)
    skew_to_vector(S_start)

    return skew_to_vector(S_start), φ_start
end


