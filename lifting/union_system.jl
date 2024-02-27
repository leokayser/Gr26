export UnionSystem, union_system

"""
   UnionSystem(F::AbstractSystem)

Documentation???
"""
struct UnionSystem <: AbstractSystem
    # The system is (f_1,...,f_r)
    F::Vector{T} where T<:AbstractSystem
    num_eq::Int64
    parameters::Vector{Variable}
end
function UnionSystem(F::Vector{T} where T<:AbstractSystem)
    num_eq = 0
    pF = parameters(F[1])

    for i in eachindex(F)
        (variable_groups(F[1]) == variable_groups(F[i]) && variables(F[1]) == variables(F[i])) || throw(
            ArgumentError(
                "Cannot create union ⋃F since the sets of variables of `F_1` and `F_$i` don't match.",
            ),
        )
        num_eq += size(F[i], 1)

        if i==1 continue end
        pFi = parameters(F[i])
        (isempty(pF) || isempty(pFi) || pF == pFi) || throw(
            ArgumentError(
                "Cannot construct a union of systems with different sets of parameters.",
            ),
        )

        pF = isempty(pFi) ? [pFi; pF] : pFi
    end
        
    UnionSystem(F, num_eq, pF)
end

Base.size(US::UnionSystem) = (US.num_eq, size(US.F[1], 2))

HomotopyContinuation.ModelKit.variables(US::UnionSystem) = variables(US.F[1])
HomotopyContinuation.ModelKit.parameters(US::UnionSystem) = US.parameters
HomotopyContinuation.ModelKit.variable_groups(US::UnionSystem) = variable_groups(US.F[1])

function Base.show(io::IO, US::UnionSystem)
    num_syss = length(US.F)
    print(io, "Union F_1∪...∪F_$num_syss:")
    for i in 1:length(US.F)
        println(io, "\nF_$i: ")
        show(io, US.F[i])
    end
end

(US::UnionSystem)(x, p = nothing) = Base.vcat([f(x, p) for f in US.F]...)
function HomotopyContinuation.ModelKit.evaluate!(
    u,
    US::UnionSystem,
    x,
    p = nothing,
)
    offset = 0
    for f in US.F
        l = size(f, 1)
        HomotopyContinuation.ModelKit.evaluate!(view(u, (offset+1):(offset+l)), f, x, p)
        offset += l
    end
end

function ModelKit.evaluate_and_jacobian!(u, U, US::UnionSystem, x, p = nothing)
    offset = 0
    for f in US.F
        l = size(f,1)
        HomotopyContinuation.ModelKit.evaluate_and_jacobian!(
            view(u, (offset[1]+1):(offset[1]+l)),
            view(U, (offset[1]+1):(offset[1]+l),:),
            f, x, p)
        offset += l
    end
end


function ModelKit.taylor!(u, v::Val{1}, US::UnionSystem, tx, p = nothing)
    offset = 0
    for f in US.F
        l = size(f,1)
        vie = view(u, (offset[1]+1):(offset[1]+l))
        tvv = TaylorVector{2}(vie)
        HomotopyContinuation.ModelKit.taylor!(
            tvv, v, f, tx, p)
        offset += l
    end
end




"""
Documentation????
"""
union_system(F::Vector{T} where T<:AbstractSystem) = UnionSystem(F)
union_system(
    F::Vector{T} where T<:System;
    compile::Union{Bool,Symbol} = HomotopyContinuation.COMPILE_DEFAULT[],
) = UnionSystem([fixed(f; compile = compile) for f in F])
