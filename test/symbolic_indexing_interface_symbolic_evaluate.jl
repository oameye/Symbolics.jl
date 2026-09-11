using Symbolics
using SymbolicIndexingInterface
using Symbolics: Differential, Operator, value

@variables t x(t) y(t)
@variables p[1:3, 1:3] q[1:3]

bar(x, p) = p * x
@register_array_symbolic bar(x::AbstractVector, p::AbstractMatrix) begin
    size = size(x)
    eltype = promote_type(eltype(x), eltype(p))
    ndims = 1
end

D = Differential(t)

expr1 = x + y + D(x)
@test isequal(symbolic_evaluate(expr1, Dict(x => 3)), 3 + y + D(3))
@test isequal(symbolic_evaluate(expr1, Dict(x => 3); operator = Operator), 3 + y + D(x))
@test isequal(symbolic_evaluate(expr1, Dict(x => 1, D(x) => 2)), y + 3)
@test value(symbolic_evaluate(expr1, Dict(x => 1, D(x) => 2, y => 3))) == 6
@test isequal(symbolic_evaluate(expr1, Dict(x => 3, y => 3x), operator = Operator), 12 + D(x))
@test value(symbolic_evaluate(expr1, Dict(x => 3, y => 3x, D(x) => 2))) == 14

expr2 = bar(q, p)
@test isequal(symbolic_evaluate(expr2, Dict(p => ones(3, 3))), bar(q, ones(3, 3)))
@test value(symbolic_evaluate(expr2, Dict(p => ones(3, 3), q => ones(3)))) == 3ones(3)

expr3 = bar(3q, 3p)
@test isequal(symbolic_evaluate(expr3, Dict(p => ones(3, 3))), bar(3q, 3ones(3, 3)))
@test value(symbolic_evaluate(expr3, Dict(p => ones(3, 3), q => ones(3)))) == 27ones(3)

expr4 = D(x) ~ 3x + y
@test isequal(symbolic_evaluate(expr4, Dict(x => 3)), D(3) ~ 9 + y)
@test isequal(symbolic_evaluate(expr4, Dict(x => 3); operator = Operator), D(x) ~ y + 9)
@test isequal(symbolic_evaluate(expr4, Dict(x => 1, D(x) => 2)), 2 ~ 3 + y)
@test isequal(symbolic_evaluate(expr4, Dict(x => 1, D(x) => 2, y => 3)), 2 ~ 6)

# General numeric symbolic scalars participate in the same indexing interface as Num.
@variables z::Number
expr5 = 1 + im * z
@test symbolic_type(typeof(z)) == ScalarSymbolic()
@test symbolic_type(typeof(expr5)) == ScalarSymbolic()
@test isequal(symbolic_evaluate(z, Dict(z => im)), im)
@test isequal(symbolic_evaluate(expr5, Dict(z => 2)), 1 + 2im)

# Regression for #1950: complex-valued matrix substitutions must remain symbolic arrays
# so variable discovery can see symbolic entries.
θ, λ = @variables θ::Real λ::Real
real_f = Symbolics.variable(:real_f; T = Symbolics.FnType{Tuple{Vararg{Number}}, Number, Nothing})(θ)
complex_f = Symbolics.variable(:complex_f; T = Symbolics.FnType{Tuple{Vararg{Number}}, Number, Nothing})(λ)
real_mat = Matrix([exp(θ) 0.0; 0.0 0.0])
complex_mat = Matrix([exp(im * λ) 0.0; 0.0 0.0])
expr6 = real_f * complex_f
real_subs = substitute(expr6, Dict(real_f => real_mat))
complex_subs = substitute(expr6, Dict(complex_f => complex_mat))
@test θ in get_variables(real_subs)
@test λ in get_variables(complex_subs)
