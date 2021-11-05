
"""
    LinearizationResult{T, X}

Represents a linear state-space system
```
ẋ = Ax + Bu
y = Cx + Du
```
obtained through a call to [`linearize`](@ref). Contains the fields `A, B, C, D, x, u, y`
"""
struct LinearizationResult{T,ET,X}
    A::T
    E::ET
    B::T
    C::T
    D::T
    x::X
    u::X
    y::X
end

function LinearizationResult(A,E,B,C,D,x,u,y)
    A,B,C,D = promote(A,B,C,D)
    x,u,y = promote(x,u,y)
    LinearizationResult(A,E,B,C,D,x,u,y)
end

"""
    linearize(
        sys::ODESystem;
        x0 = nothing,
        u = unbound_inputs(sys),
        y = unbound_outputs(sys),
        numeric = false,
    )

Linaerize `sys` around operating point `x0`. `x0` can be `nothing` for a symbolic linearization, or a `Dict` that maps states and parameters to values.

Returns an instance of [`LinearizationResult`](@ref)
representing a linear state-space system
```
ẋ = Ax + Bu
y = Cx + Du
```

If `numeric = true`, symbolic parameters are replaced with their values specified in `x0`.
Parameters that are not specified in `x0` will use defaults if available, otherwise an error is thrown.
"""
function linearize_symbolic(
    sys::ODESystem;
    x0 = nothing,
    u = unbound_inputs(sys),
    y = unbound_outputs(sys),
    numeric = false,
)
    sysr = inconsistent_simplification(sys; simplify = false)
    x = ModelingToolkit.collect_differential_variables(sysr) |> collect
    dynamics = [e.rhs for e in equations(sysr) if isdifferential(e.lhs)]
    oexprs = output_expressions(sysr, y)
    A = Symbolics.jacobian(dynamics, x)
    B = Symbolics.jacobian(dynamics, u)
    C = Symbolics.jacobian(oexprs, x)
    D = Symbolics.jacobian(oexprs, u)

    @assert size(A, 1) == size(A, 2) == length(x)
    @assert size(A, 1) == size(B, 1)
    @assert size(B, 2) == length(u)
    @assert size(C) == (length(y), length(x))
    @assert size(D) == (length(y), length(u))

    linres = LinearizationResult(A, I(size(A,1)), B, C, D, x, u, y)

    if numeric || x0 !== nothing
        sym2val = ModelingToolkit.defaults(sysr)
        if x0 isa AbstractDict
            sym2val = merge(sym2val, x0) # NOTE: it's important to have x0 last for x0 to take precedence
        end
        linres = substitute(linres, sym2val)
    end

    # x = states(sysr, states(sysr)) # redo this to get the namespace right for later use
    # u = filter(isinput, x)
    # x = filter(!isinput, x)
    LinearizationResult(A, I(size(A,1)), B, C, D, x, u, y)
end

function linearize(
    sys::ODESystem;
    x0 = nothing,
    u = unbound_inputs(sys),
    y = unbound_outputs(sys),
)
    sym2val = ModelingToolkit.defaults(sys)
    if x0 isa AbstractDict
        sym2val = merge(sym2val, x0) # NOTE: it's important to have x0 last for x0 to take precedence
    end

    A,E,B = dss_fd(sys; x0 = sym2val, p = sym2val)

    @assert size(A, 1) == size(A, 2) == length(x)
    @assert size(A, 1) == size(B, 1)
    @assert size(B, 2) == length(u)
    @assert size(C) == (length(y), length(x))
    @assert size(D) == (length(y), length(u))

    linres = LinearizationResult(A, E, B, C, D, x, u, y)
end

function get_AEB(sysr, AB, CD = nothing)
    x = states(sysr)
    size(AB, 2) == length(x) || throw(ArgumentError("Width of AE does not match length(states(sys))"))
    dx = collect_differential_variables(sysr)

    setdx = Set(dx)
    diffeq_inds = map(eq-> ModelingToolkit.isdifferential(eq.lhs), equations(sysr))
    state_inds = map(x->x ∈ setdx, x)
    control_inds = map(x->isinput(x) && !is_bound(sysr, x), x)
    @assert !any(state_inds .& control_inds) # cannot be both state and control

    nx = length(x) - count(control_inds)
    A = AB[:, .!control_inds]
    @assert size(A) == (nx, nx)
    B = AB[:, control_inds]
    E = zeros(size(A,1), size(A,1))
    for i = findall(diffeq_inds)
        E[i,i] = 1
    end
    if CD === nothing
        (; A, E, B)
    else
        C = CD[:, .!control_inds]
        D = CD[:, control_inds]
        return (; A, E, B, C, D)
    end
end

using ForwardDiff


function dss_fd(sys::ODESystem;
    x0 = nothing,
    p = nothing,
    )
    sysr = inconsistent_simplification(sys; simplify = true)
    if x0 === nothing
        x0 = ModelingToolkit.varmap_to_vars(ModelingToolkit.defaults(sysr), states(sysr))
    end
    if p === nothing
        p = ModelingToolkit.varmap_to_vars(ModelingToolkit.defaults(sysr), parameters(sysr))
    end
    fun = ODEFunction{true}(sysr)
    AB = ForwardDiff.jacobian(x->fun(x, p, 1), x0)
    out = zeros(size(x0))
    ny = length(observed(sysr))
    out = zeros(ny)
    CD = ForwardDiff.jacobian(x->fun.observed(out, x, p, 1), x0)
    get_AEB(sysr, AB, CD)
end

function dss(sys::ODESystem)
    sysr = inconsistent_simplification(sys; simplify = true)
    sym2val = ModelingToolkit.defaults(sysr)
    jac0 = calculate_jacobian(sysr)
    AB = substitute.(jac0, Ref(sym2val)) .|> Float64
    get_AEB(sysr, AB)
end

"""
    inconsistent_simplification(sys; simplify = false)

Perform [`structural_simplify`](@ref) without consistency check.
"""
function inconsistent_simplification(sys; simplify = false)
    sysr = initialize_system_structure(alias_elimination(sys))
    if sysr isa ODESystem
        sysr = dae_index_lowering(sysr)
    end
    sysr = tearing(sysr, simplify = simplify)
end

"""
    substitute(linres, sym2val::Dict, eltype=Float64)

Substitute symbolic variables in the system representation `A,B,C,D` with numeric values.
`sym2val` is a dict mapping symbolic variables to `Number`s.
"""
function SymbolicUtils.substitute(
    linres::LinearizationResult,
    sym2val::AbstractDict,
    ::Type{FT} = Float64,
) where {FT}
    function tofloat(x)
        xv = value(x)
        xv isa Real || error(
            "There were symbolic variables remaining after substitution. This means that `sym2val` did not contain a map for all relevant values. If you called linearize, make sure that there are either defaults for all relevant parameters and states, or that those parameters and states that lack defaults are provided in the argument `x0`. The problematic expression is $(xv)",
        )
        FT(xv)
    end
    @unpack A, B, C, D = linres
    A = substitute.(A, Ref(sym2val)) .|> tofloat
    B = substitute.(B, Ref(sym2val)) .|> tofloat
    C = substitute.(C, Ref(sym2val)) .|> tofloat
    D = substitute.(D, Ref(sym2val)) .|> tofloat
    LinearizationResult(
        A,
        B,
        C,
        D,
        ntuple(i -> getfield(linres, 4 + i), fieldcount(typeof(linres)) - 4)...,
    )
end

# """
#     differential_states(sys)

# Return a vector of differential states, i.e., states `xᵢ` that are associated with a differential equation
# `dxᵢ / dt = f(x, u)`
# """
# function differential_states(sys)
#     diffvars = Set([e.lhs.arguments[] for e in equations(sys) if isdifferential(e.lhs)])
#     x = [map(eq -> eq.lhs, observed(sys)); states(sys)]
#     filter(x -> x ∈ diffvars, x)
# end

"""
    output_expressions(sysr, ys::Vector)

Return the lhs of all equations defining outputs.
"""
function output_expressions(sysr, ys::AbstractVector)
    sety = Set(ys)
    oexprs = Num[]
    eqs = equations(sysr)
    oeqs = observed(sysr)
    aeqs = [oeqs; eqs]
    for eq in aeqs
        vars  = get_equationvars(eq)
        numys = count(v -> v ∈ sety, vars)
        numys == 0 && continue
        numys > 1 && error("Found more than one linearization output in the same equation.")
        y  = vars[findfirst(v -> v ∈ sety, vars)] # The particular output variable in this eq
        ex = Symbolics.solve_for(eq, y) # get the expression equalling y
        push!(oexprs, ex)
    end
    # pass through expressions and substitute observed variables with non-observed variables until no observed variables remain

    oexprs
end

