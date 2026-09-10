# Experimental replacement for the legacy Num/Complex methods in num.jl. Ordinary
# arithmetic with a numerical complex value should stay in the raw symbolic algebra and
# select its wrapper from the result symtype; it should not Cartesian-expand into
# Complex{Num}. The final change should delete the old methods rather than redefine them.
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
