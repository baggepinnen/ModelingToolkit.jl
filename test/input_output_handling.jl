using ModelingToolkit, Symbolics, Test
using ModelingToolkit: get_namespace, has_var, inputs, outputs, is_bound, bound_inputs, unbound_inputs

@parameters tv
D = Differential(tv)
@variables x(tv) u(tv) [input=true]
@named sys = ODESystem([D(x) ~ -x + u], tv) # both u and x are unbound
@named sys2 = ODESystem([D(x) ~ -sys.x], tv, systems=[sys]) # this binds sys.x in the context of sys2, sys2.x is still unbound
@named sys3 = ODESystem([D(x) ~ -sys.x + sys.u], tv, systems=[sys]) # This binds both sys.x and sys.u
@test has_var(x ~ 1, x)
@test has_var(1 ~ x, x)
@test has_var(x + x, x)
@test !has_var(2 ~ 1, x)

@test get_namespace(x) == ""
@test get_namespace(sys.x) == "sys"
@test get_namespace(sys2.x) == "sys2"
@test get_namespace(sys2.sys.x) == "sys2₊sys"

@test !is_bound(sys, u)
@test !is_bound(sys, x)
@test !is_bound(sys, sys.u)
@test  is_bound(sys2, sys.x)
@test !is_bound(sys2, sys.u)
@test !is_bound(sys2, sys2.sys.u)

@test is_bound(sys3, sys.u) # I would like to write sys3.sys.u here but that's not how the variable is stored in the equations
@test is_bound(sys3, sys.x)

@test isequal(inputs(sys), [u])
@test isequal(inputs(sys2), [sys.u])

@test isempty(bound_inputs(sys))
@test isequal(unbound_inputs(sys), [u])

@test isempty(bound_inputs(sys2))
@test isequal(unbound_inputs(sys2), [sys.u])

@test isequal(bound_inputs(sys3), [sys.u])
@test isempty(unbound_inputs(sys3))