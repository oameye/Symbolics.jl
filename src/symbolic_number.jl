# Atomic wrapper for symbolic scalar expressions whose symbolic type is numeric but not
# necessarily real. Unlike `Complex{Num}`, this wrapper does not decompose an expression
# into real and imaginary components. The wrapped `BasicSymbolic` remains one expression
# tree.
@symbolic_wrap struct SymbolicNumber <: Number
    val::BasicSymbolic{VartypeT}

    function SymbolicNumber(ex::BasicSymbolic{VartypeT})
        @assert symtype(ex) <: Number
        return new(Const{VartypeT}(ex))
    end

    function SymbolicNumber(ex::Number)
        return new(Const{VartypeT}(unwrap(ex)))
    end
end

SymbolicNumber(x::SymbolicNumber) = x

SymbolicUtils.unwrap(x::SymbolicNumber) = x.val
SU.infer_vartype(::Type{SymbolicNumber}) = VartypeT
SymbolicUtils.symtype(x::SymbolicNumber) = symtype(unwrap(x))

# Route symbolic arithmetic through the raw BasicSymbolic algebra and only choose the
# wrapper after `promote_symtype` has determined the mathematical result domain. `//` is
# deliberately excluded: it constructs an exact Rational and is not generic division.
SymbolicUtils.@number_methods(
    SymbolicNumber,
    wrap(f(unwrap(a))),
    wrap(f(unwrap(a), unwrap(b))),
    [conj, real, imag, transpose, //],
)

Base.conj(x::SymbolicNumber) = wrap(conj(unwrap(x)))
Base.real(x::SymbolicNumber) = wrap(real(unwrap(x)))
Base.imag(x::SymbolicNumber) = wrap(imag(unwrap(x)))
Base.transpose(x::SymbolicNumber) = wrap(transpose(unwrap(x)))
Base.adjoint(x::SymbolicNumber) = wrap(adjoint(unwrap(x)))
# `@number_methods` defines `^(::SymbolicNumber, ::Real)`, which intersects Base's
# integer/rational power methods. Keep those powers on the symbolic algebra explicitly.
Base.:^(x::SymbolicNumber, p::Integer) = wrap(unwrap(x)^p)
Base.:^(x::SymbolicNumber, p::Rational) = wrap(unwrap(x)^p)
# Base has a dedicated `ℯ ^ ::Number` method which intersects the generic symbolic
# exponent methods emitted above. Preserve the canonical exponential representation.
Base.:^(::Irrational{:ℯ}, x::SymbolicNumber) = wrap(exp(unwrap(x)))

# `polygamma(::Integer, ::Number)` in SpecialFunctions intersects the generic symbolic
# binary-function methods. This exact intersection keeps integer orders on the symbolic
# expression path without broadening the dispatch surface.
SpecialFunctions.polygamma(m::Integer, x::SymbolicNumber) =
    wrap(SpecialFunctions.polygamma(m, unwrap(x)))

# Base implements `cis(::Real)` through `sincos` followed by explicit `Complex`
# construction. That is appropriate for numerical values but would reintroduce Cartesian
# storage for `Num`. Canonically lower symbolic `cis` to the equivalent atomic scalar
# expression instead.
Base.cis(x::Num) = wrap(exp(im * unwrap(x)))

Base.iszero(x::SymbolicNumber) = SymbolicUtils._iszero(unwrap(x))
Base.isone(x::SymbolicNumber) = SymbolicUtils._isone(unwrap(x))
Base.zero(::SymbolicNumber) = SymbolicNumber(0)
Base.zero(::Type{SymbolicNumber}) = SymbolicNumber(0)
Base.one(::SymbolicNumber) = SymbolicNumber(1)
Base.one(::Type{SymbolicNumber}) = SymbolicNumber(1)

# `SymbolicNumber` is the wide numeric wrapper. Ordinary real values mixed with `Num`
# continue to promote to `Num` via `num.jl`; only the explicit Num/complex edge in
# `complex.jl` widens a real symbolic value to `SymbolicNumber`.
Base.promote_rule(::Type{SymbolicNumber}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: Number} = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{T}) where {T <: Number} = SymbolicNumber
# Exact intersections with Base promotion rules keep Aqua ambiguity-free.
Base.promote_rule(::Type{Bool}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: AbstractIrrational} =
    SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{SymbolicNumber}) = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{Num}) = SymbolicNumber
Base.convert(::Type{SymbolicNumber}, x::Number) = SymbolicNumber(x)

# Wrappers are representation boundaries, not distinct symbolic identities. Matching the
# wrapped expression's hash and `isequal` semantics is required by generic substitution,
# which recursively visits raw `BasicSymbolic` nodes while users naturally provide wrapped
# variables as dictionary keys.
Base.hash(x::SymbolicNumber, h::UInt) = hash(unwrap(x), h)::UInt
Base.isequal(a::SymbolicNumber, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::SymbolicNumber, b::BasicSymbolic) = isequal(unwrap(a), b)
Base.isequal(a::BasicSymbolic, b::SymbolicNumber) = isequal(a, unwrap(b))
Base.isequal(a::SymbolicNumber, b::Num) = isequal(unwrap(a), unwrap(b))
Base.isequal(a::Num, b::SymbolicNumber) = isequal(unwrap(a), unwrap(b))

function Base.show(io::IO, x::SymbolicNumber)
    warn_load_latexify()
    show(io, unwrap_const(unwrap(x)))
end

# Generic symbolic utilities must select the wrapper from the transformed expression's
# resulting symtype. In particular, a simplification is allowed to narrow a
# complex-capable expression to a provably real `Num`.
SymbolicUtils.simplify(x::SymbolicNumber; kw...) = wrap(SymbolicUtils.simplify(unwrap(x); kw...))
SymbolicUtils.simplify_fractions(x::SymbolicNumber; kw...) = wrap(SymbolicUtils.simplify_fractions(unwrap(x); kw...))
SymbolicUtils.expand(x::SymbolicNumber) = wrap(SymbolicUtils.expand(unwrap(x)))
SymbolicUtils.Code.toexpr(x::SymbolicNumber) = SymbolicUtils.Code.toexpr(unwrap(x))
SymbolicUtils.setmetadata(x::SymbolicNumber, t, v) = wrap(SymbolicUtils.setmetadata(unwrap(x), t, v))
SymbolicUtils.getmetadata(x::SymbolicNumber, t) = SymbolicUtils.getmetadata(unwrap(x), t)
SymbolicUtils.hasmetadata(x::SymbolicNumber, t) = SymbolicUtils.hasmetadata(unwrap(x), t)
Broadcast.broadcastable(x::SymbolicNumber) = x
SymbolicUtils.scalarize(x::SymbolicNumber) = wrap(SymbolicUtils.scalarize(unwrap(x)))

function SymbolicUtils.search_variables!(buffer, expr::SymbolicNumber; kw...)
    SymbolicUtils.search_variables!(buffer, unwrap(expr); kw...)
end

SymbolicIndexingInterface.symbolic_type(::Type{SymbolicNumber}) = ScalarSymbolic()
SymbolicIndexingInterface.hasname(x::SymbolicNumber) = hasname(unwrap(x))
SymbolicIndexingInterface.getname(x::SymbolicNumber) = getname(unwrap(x))
function SymbolicIndexingInterface.symbolic_evaluate(x::SymbolicNumber, d::Dict; kw...)
    SymbolicIndexingInterface.symbolic_evaluate(unwrap(x), d; kw...)
end

function (s::SymbolicUtils.Substituter)(x::SymbolicNumber)
    wrap(s(unwrap(x)))
end
