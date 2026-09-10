include("symbolic_number.jl")

SymbolicUtils.promote_symtype(::typeof(imag), ::Type{Complex{T}}) where {T} = T
Base.promote_rule(::Type{Complex{T}}, ::Type{Num}) where {T <: Real} = SymbolicNumber
Base.promote_rule(::Type{Num}, ::Type{Complex{T}}) where {T <: Real} = SymbolicNumber

# A numerical complex coefficient does not imply Cartesian symbolic storage. Build the
# operation in BasicSymbolic and select Num/SymbolicNumber from its resulting symtype.
for C in (Complex, Complex{Bool})
    @eval begin
        Base.:+(x::Num, z::$C) = wrap(unwrap(x) + z)
        Base.:+(z::$C, x::Num) = wrap(z + unwrap(x))
        Base.:-(x::Num, z::$C) = wrap(unwrap(x) - z)
        Base.:-(z::$C, x::Num) = wrap(z - unwrap(x))
        Base.:*(x::Num, z::$C) = wrap(unwrap(x) * z)
        Base.:*(z::$C, x::Num) = wrap(z * unwrap(x))
        Base.:/(x::Num, z::$C) = wrap(unwrap(x) / z)
        Base.:/(z::$C, x::Num) = wrap(z / unwrap(x))
    end
end

# `Complex{Num}` remains a supported explicit Cartesian representation, but it is no
# longer selected by `wrap` for symbolic complex scalars. `SymbolicNumber <: Number`
# handles those atomically via the generic symbolic-wrapper dispatch.
is_wrapper_type(::Type{Complex{Num}}) = true
wraps_type(::Type{Complex{Num}}) = Complex{Real}
iswrapped(::Complex{Num}) = true

function SymbolicUtils.unwrap(a::Complex{<:Num})
    re, img = unwrap(real(a)), unwrap(imag(a))
    if SymbolicUtils.isconst(re) && SymbolicUtils.isconst(img)
        return Const{VartypeT}(complex(unwrap_const(re), unwrap_const(img)))
    end
    if iscall(re) && operation(re) === real && iscall(img) && operation(img) === imag && isequal(arguments(re)[1], arguments(img)[1])
        return arguments(re)[1]
    end
    sT = promote_type(symtype(re), symtype(img))
    return Term{VartypeT}(complex, SymbolicUtils.ArgsT{vartype(re)}((re, img)); type = Complex{sT}, shape = SymbolicUtils.ShapeVecT())
end

SymbolicUtils.infer_vartype(::Type{Complex{Num}}) = VartypeT

function Base.Complex{Num}(x::BasicSymbolic{VartypeT})
    Complex{Num}(wrap(real(x)), wrap(imag(x)))
end

function Base.show(io::IO, a::Complex{Num})
    rr = unwrap(real(a))
    ii = unwrap(imag(a))

    if iscall(rr) && (operation(rr) === real) &&
        iscall(ii) && (operation(ii) === imag) &&
        isequal(arguments(rr)[1], arguments(ii)[1])

        return print(io, arguments(rr)[1])
    end

    show(io, real(a) + im * imag(a))
end

function (s::SymbolicUtils.Substituter)(x::Complex{Num})
    Complex{Num}(s(real(x)), s(imag(x)))
end
