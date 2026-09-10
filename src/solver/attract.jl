"""
    remove_var_factor(lhs, var)

Divides out factors containing the variable from both sides of an equation expression.
"""
function remove_var_factor(lhs, var)
    lhs = unwrap(lhs)
    if ismul(lhs)
        new_args = [arg for arg in arguments(lhs) if !contains_var(arg, var)]
        isempty(new_args) && return one(lhs)
        return prod(new_args)
    end
    return lhs
end

# NOTE: file content before attract_exponential is intentionally preserved semantically by
# the branch; the relevant complex fix is in the rule below.

"""
    attract_exponential(lhs, var)

Rewrites ``a*b^f(x) + c*d^g(x)`` into
``f(x) * log(b) - g(x) * log(d) + log(-a/c)``.
"""
function attract_exponential(lhs, var)
    lhs = unwrap(lhs)
    contains_var(arg) = n_occurrences(arg, var) > 0

    r_addexpon = Vector{Any}()

    #! format: off
    push!(r_addexpon, @acrule (~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -1))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, -~a))
    push!(r_addexpon, @acrule (~a)*(~b)^(~f::(contains_var)) + (~c)*(~d)^(~g::(contains_var)) => ~f*term(slog, ~b) - ~g*term(slog, ~d) + term(slog, sdiv(-(~a), ~c)))
    #! format: on

    lhs = expand(simplify(
        lhs, rewriter = SymbolicUtils.Postwalk(SymbolicUtils.Chain(r_addexpon))))

    return expand(lhs)
end
