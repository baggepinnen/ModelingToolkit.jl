using Symbolics: get_variables
"""
    inputs(sys)

Return all variables that mare marked as inputs. See also [`unbound_inputs`](@ref)
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref)
"""
inputs(sys) = filter(isinput, states(sys))

"""
    outputs(sys)

Return all variables that mare marked as outputs. See also [`unbound_outputs`](@ref)
See also [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function outputs(sys)
    o = observed(sys)
    rhss = [eq.rhs for eq in o]
    lhss = [eq.lhs for eq in o]
    unique([
        filter(isoutput, states(sys))
        filter(x -> x isa Term && isoutput(x), rhss) # observed can return equations with complicated expressions, we are only looking for single Terms
        filter(x -> x isa Term && isoutput(x), lhss)
    ])
end

"""
    bound_inputs(sys)

Return inputs that are bound within the system, i.e., internal inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_inputs(sys) = filter(x->is_bound(sys, x), inputs(sys))

"""
    unbound_inputs(sys)

Return inputs that are not bound within the system, i.e., external inputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_inputs(sys) = filter(x->!is_bound(sys, x), inputs(sys))

"""
    bound_outputs(sys)

Return outputs that are bound within the system, i.e., internal outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
bound_outputs(sys) = filter(x->is_bound(sys, x), outputs(sys))

"""
    unbound_outputs(sys)

Return outputs that are not bound within the system, i.e., external outputs
See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
unbound_outputs(sys) = filter(x->!is_bound(sys, x), outputs(sys))

"""
    is_bound(sys, u)

Determine whether or not input/output variable `u` is "bound" within the system, i.e., if it's to be considered internal to `sys`.
A variable/signal is considered bound if it appears in an equation together with variables from other subsystems.
The typical usecase for this function is to determine whether the input to an IO component is connected to another component,
or if it remains an external input that the user has to supply before simulating the system. 

See also [`bound_inputs`](@ref), [`unbound_inputs`](@ref), [`bound_outputs`](@ref), [`unbound_outputs`](@ref)
"""
function is_bound(sys, u, stack=[])
    #=
    For observed quantities, we check if a variable is connected to something that is bound to something further out. 
    In the following scenario
    julia> observed(syss)
        2-element Vector{Equation}:
        sys₊y(tv) ~ sys₊x(tv)
        y(tv) ~ sys₊x(tv)
    sys₊y(t) is bound to the outer y(t) through the variable sys₊x(t) and should thus return is_bound(sys₊y(t)) = true.
    When asking is_bound(sys₊y(t)), we know that we are looking through observed equations and can thus ask 
    if var is bound, if it is, then sys₊y(t) is also bound. This can lead to an infinite recursion, so we maintain a stack of variables we have previously asked about to be able to break cycles
    =#
    u ∈ Set(stack) && return false # Cycle detected
    eqs = equations(sys)
    eqs = filter(eq->has_var(eq, u), eqs) # Only look at equations that contain u
    # isout = isoutput(u)
    for eq in eqs
        vars = get_equationvars(eq)
        for var in vars
            var === u && continue
            if !same_or_inner_namespace(u, var)
                return true
            end
        end
    end
    # Look through observed equations as well
    oeqs = observed(sys)
    oeqs = filter(eq->has_var(eq, u), oeqs) # Only look at equations that contain u
    for eq in oeqs
        vars = get_equationvars(eq)
        for var in vars
            var === u && continue
            if !same_or_inner_namespace(u, var)
                return true
            end
            if is_bound(sys, var, [stack; u]) && !inner_namespace(u, var) # The variable we are comparing to can not come from an inner namespace, binding only counts outwards
                return true
            end
        end
    end
    false
end

"""
    same_or_inner_namespace(u, var)

Determine whether or not `var` is in the same namespace as `u`, or a namespace internal to the namespace of `u`.
Example: `sys.u ~ sys.inner.u` will bind `sys.inner.u`, but `sys.u` remains an unbound, external signal. The namepsaced signal `sys.inner.u` lives in a namspace internal to `sys`.
"""
function same_or_inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv ||           # namespaces are the same
        startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namepsace to nu
        occursin('₊', var) && !occursin('₊', u) # or u is top level but var is internal
end

function inner_namespace(u, var)
    nu = get_namespace(u)
    nv = get_namespace(var)
    nu == nv && return false
    startswith(nv, nu) || # or nv starts with nu, i.e., nv is an inner namepsace to nu
        occursin('₊', var) && !occursin('₊', u) # or u is top level but var is internal
end

"""
    get_namespace(x)

Return the namespace of a variable as a string. If the variable is not namespaced, the string is empty.
"""
function get_namespace(x)
    sname = string(x)
    parts = split(sname, '₊')
    if length(parts) == 1
        return ""
    end
    join(parts[1:end-1], '₊')
end

"""
    has_var(eq, x)

Determine whether or not an equation or expression contains variable `x`.
"""
function has_var(eq::Equation, x)
    has_var(eq.rhs, x) || has_var(eq.lhs, x)
end

has_var(ex, x) = x ∈ Set(get_variables(ex))

get_equationvars(eq) = [get_variables(eq.rhs); get_variables(eq.lhs)]

"""
    linearize(sys::ODESystem; numeric=false, x0=nothing)

Linaerize `sys` around operating point `x0`. `x0` can be `nothing` for a symbolic linearization, or a `Dict` that maps states and parameters to values.

Returns a named tuple with fields `A,B,x,u`
representing a linear state-space system
```
ẋ = Ax + Bu
```

If `numeric = true`, symbolics parameters are replaced with their default values and the matrices `A,B` will have element type `Float64`, only pass this option if all parameters and states have default values or have values in `x0`.
"""
function linearize(sys::ODESystem; x0=nothing, u = unbound_inputs(sys), y = unbound_outputs(sys), numeric=false)
    sysr = inconsistent_simplification(sys; simplify=false)
    x = differential_states(sysr)
    # Set(x) == Set(states(sysr)) || error("Algebraic variables remaining, this is currently not supported in linearize.")
    dynamics = [e.rhs for e in equations(sysr) if isdifferential(e.lhs)]
    oexprs = output_expressions(sysr, y)
    A = Symbolics.jacobian(dynamics, x)
    B = Symbolics.jacobian(dynamics, u) 
    C = Symbolics.jacobian(oexprs, x) 
    D = Symbolics.jacobian(oexprs, u) 

    if numeric || x0 !== nothing
        sym2val = ModelingToolkit.defaults(sysr)
        if x0 isa AbstractDict
            sym2val = merge(sym2val, x0) # NOTE: it's important to have x0 last for x0 to take precedence
        end
        A,B,C,D = numeric_linearization(A, B, C, D, sym2val)
    end
    @assert size(A,1) == size(A,2) == length(x)
    @assert size(A,1) == size(B,1)
    @assert size(B,2) == length(u)
    @assert size(C) == (length(y), length(x))
    @assert size(D) == (length(y), length(u))
    # x = states(sysr, states(sysr)) # redo this to get the namespace right for later use
    # u = filter(isinput, x)
    # x = filter(!isinput, x)
    (; A, B, C, D, x, u)
end

function inconsistent_simplification(sys; simplify=false)
    sysr = initialize_system_structure(alias_elimination(sys))
    if sysr isa ODESystem
        sys = dae_index_lowering(sysr)
    end
    sysr = tearing(sys, simplify=simplify)
end

function numeric_linearization(A, B, C, D, sym2val)
    function tofloat(x)
        xv = value(x)
        xv isa Real || error("There were symbolic variables remaining after substitution. This means that `sym2val` did not contain a map for all relevant values. If you called linearize, make sure that there are either defaults for all relevant parameters and states, or that those parameters and states that lack defaults are provided in the argument `x0`. The problematic expression is $(xv)")
        Float64(xv)
    end
    A = substitute.(A, Ref(sym2val)) .|> tofloat
    B = substitute.(B, Ref(sym2val)) .|> tofloat
    C = substitute.(C, Ref(sym2val)) .|> tofloat
    D = substitute.(D, Ref(sym2val)) .|> tofloat
    (; A, B, C, D)
end

function differential_states(sys)
    diffvars = Set([e.lhs.arguments[] for e in equations(sys) if isdifferential(e.lhs)])
    x = [map(eq->eq.lhs, observed(sys)); states(sys)]

    filter(x->x ∈ diffvars, x)
end

function output_expressions(sysr, ys)
    sety = Set(ys)
    oexprs = Num[]
    eqs = [observed(sysr); equations(sysr)]
    for eq in eqs
        vars = get_equationvars(eq)
        numys = count(v->v ∈ sety, vars)
        numys == 0 && continue
        numys > 1 && error("Found more than one linearization output in the same equation.")
        y = vars[findfirst(v->v ∈ sety, vars)] # The particular output variable in this eq
        ex = Symbolics.solve_for(eq, y) # get the expression equalling y
        push!(oexprs, ex)
    end
    oexprs
end