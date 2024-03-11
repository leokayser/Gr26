

using DelimitedFiles
Γ = randomSApoints(14)

using Pkg
Pkg.add("JLD")
using JLD
JLD.save("data.jld", "data", Γ)
JLD.load("data.jld")["data"] == Γ


function normalize_SApoints(Γ)
    ind = max_independent_points(Γ)
    Γ = hcat(Γ[:,ind],Γ[:,deleteat!([i for i=1:14],ind)])
    λ = self_duality_control(Γ)
    λ_normalizing = diagm(vcat([sqrt(λ[i,i]) for i=1:7],[sqrt(-λ[i,i]) for i=8:14]))
    Γ_scaled = Γ * λ_normalizing
    O = inv(Γ_scaled[:,1:7]) * Γ_scaled[:,8:14] 
    Γnorm = hcat( diagm([1 for i=1:7]) , O )
    return (Γnorm, Γ_scaled[:,1:7], λ_normalizing)
end