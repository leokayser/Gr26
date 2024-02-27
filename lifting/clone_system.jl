export CloneSystem, clone_sys

struct CloneSystem{S1<:AbstractSystem} <: AbstractSystem
    f::S1
end
#function CloneSystem(f::AbstractSystem)
#    CloneSystem(f)
#end

Base.size(C::CloneSystem) = size(C.f)

ModelKit.variables(F::CloneSystem) = variables(F.f)
ModelKit.parameters(F::CloneSystem) = parameters(F.f)
ModelKit.variable_groups(F::CloneSystem) = variable_groups(F.f)

function Base.show(io::IO, C::CloneSystem)
    println(io, "Clone of:")
    show(io, C.f)
end

(C::CloneSystem)(x, p = nothing) = C.f(x, p)
function ModelKit.evaluate!(
    u,
    C::CloneSystem,
    x,
    p = nothing,
)
    evaluate!(u, C.f, x, p)
end

function ModelKit.evaluate_and_jacobian!(u, U, C::CloneSystem, x, p = nothing)
    evaluate_and_jacobian!(u, U, C.f, x, p)
    nothing
end

function ModelKit.taylor!(u, v::Val, C::CompositionSystem, tx, p = nothing)
    taylor!(u, v, C.f, tx, p)
end

clone_sys(
    f::AbstractSystem;
    compile::Union{Bool,Symbol} = COMPILE_DEFAULT[],
) = CloneSystem(f)
clone_sys(f::System; compile::Union{Bool,Symbol} = HomotopyContinuation.COMPILE_DEFAULT[]) =
    CloneSystem(fixed(f; compile = compile))