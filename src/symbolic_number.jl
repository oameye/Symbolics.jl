export SymbolicNumber

"""
    SymbolicNumber

Atomic wrapper for symbolic scalar expressions whose symbolic type is numeric but not
necessarily real. Unlike `Complex{Num}`, this wrapper does not decompose an expression into
real and imaginary components. The wrapped `BasicSymbolic` remains a single expression tree.
"""
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
# wrapper after `promote_symtype` has determined the mathematical result domain.
SymbolicUtils.@number_methods(
    SymbolicNumber,
    wrap(f(unwrap(a))),
    wrap(f(unwrap(a), unwrap(b))),
    [conj, real, imag, transpose],
)

Base.conj(x::SymbolicNumber) = wrap(conj(unwrap(x)))
Base.real(x::SymbolicNumber) = wrap(real(unwrap(x)))
Base.imag(x::SymbolicNumber) = wrap(imag(unwrap(x)))
Base.transpose(x::SymbolicNumber) = wrap(transpose(unwrap(x)))
Base.adjoint(x::SymbolicNumber) = wrap(adjoint(unwrap(x)))

# Unparameterized `Complex(re, im)` is ordinary scalar construction. When both components
# are symbolic real scalars, keep the result as one symbolic scalar instead of eagerly
# changing representation to `Complex{Num}`. Users that explicitly need the Cartesian
# representation can still request `Complex{Num}(re, im)`.
function Base.Complex(re::Num, img::Num)
    return wrap(unwrap(re) + im * unwrap(img))
end

Base.iszero(x::SymbolicNumber) = SymbolicUtils._iszero(unwrap(x))
Base.isone(x::SymbolicNumber) = SymbolicUtils._isone(unwrap(x))
Base.zero(::SymbolicNumber) = SymbolicNumber(0)
Base.zero(::Type{SymbolicNumber}) = SymbolicNumber(0)
Base.one(::SymbolicNumber) = SymbolicNumber(1)
Base.one(::Type{SymbolicNumber}) = SymbolicNumber(1)

Base.promote_rule(::Type{T}, ::Type{SymbolicNumber}) where {T <: Number} = SymbolicNumber
Base.promote_rule(::Type{SymbolicNumber}, ::Type{T}) where {T <: Number} = SymbolicNumber
Base.convert(::Type{SymbolicNumber}, x::Number) = SymbolicNumber(x)

Base.hash(x::SymbolicNumber, h::UInt) = hash(unwrap(x), h)::UInt

function Base.show(io::IO, x::SymbolicNumber)
    warn_load_latexify()
    show(io, unwrap_const(unwrap(x)))
end

function SymbolicUtils.search_variables!(buffer, expr::SymbolicNumber; kw...)
    SymbolicUtils.search_variables!(buffer, unwrap(expr); kw...)
end

SymbolicIndexingInterface.symbolic_type(::Type{SymbolicNumber}) = ScalarSymbolic()
SymbolicIndexingInterface.hasname(x::SymbolicNumber) = hasname(unwrap(x))
SymbolicIndexingInterface.getname(x::SymbolicNumber) = getname(unwrap(x))

function (s::SymbolicUtils.Substituter)(x::SymbolicNumber)
    wrap(s(unwrap(x)))
end
