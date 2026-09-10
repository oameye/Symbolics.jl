using Symbolics, Test
using SymbolicUtils: metadata, symtype, unwrap_const
using Symbolics: unwrap
using SymbolicIndexingInterface: getname, hasname

@variables a b::Real z::Complex (Z::Complex)[1:10]

@testset "atomic complex scalar types" begin
    @test a isa Num
    @test b isa Num
    @test z isa Symbolics.SymbolicNumber
    @test Z[1] isa Symbolics.SymbolicNumber

    for x in (z, Z[1], z + a, z * a, z^2, z / z)
        @test x isa Symbolics.SymbolicNumber
        @test symtype(unwrap(x)) <: Number
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa Symbolics.SymbolicNumber
    end

    # Explicit Cartesian representation remains available as a compatibility path.
    cart = Complex{Num}(a, b)
    @test cart isa Complex{Num}
    @test real(cart) === a
    @test imag(cart) === b
end

@testset "literal imaginary unit stays literal" begin
    ai = a * im
    ia = im * a

    @test ai isa Symbolics.SymbolicNumber
    @test ia isa Symbolics.SymbolicNumber
    @test symtype(unwrap(ai)) <: Number
    @test symtype(unwrap(ia)) <: Number
    @test !isdefined(Symbolics, :IM)

    @test repr(ai) == "a*im" || repr(ai) == "im*a"
    @test repr(ia) == "a*im" || repr(ia) == "im*a"

    @test Symbolics.value(substitute(ai, Dict(a => 2.0))) == 2.0im
    @test Symbolics.value(substitute(ia, Dict(a => 2.0))) == 2.0im
end

@testset "substitution" begin
    # issue #314, but without changing representation to Complex{Num}
    bi = a + a * im
    @test bi isa Symbolics.SymbolicNumber

    bs = substitute(bi, Dict(a => 1.0))
    @test bs isa Symbolics.SymbolicNumber
    bv = unwrap_const(Symbolics.value(bs))
    @test bv == 1.0 + 1.0im
    @test typeof(bv) == ComplexF64
end

@testset "compact representation" begin
    @test repr(z) == "z"
    @test repr(a + b * im) == "a + b*im"
    @test repr(exp(im * a)) == "exp(im*a)" || repr(exp(im * a)) == "exp(a*im)"
end

@testset "metadata" begin
    # Complex variables are now single symbolic variables, so their metadata lives on the
    # atomic tree rather than being duplicated onto synthetic `.re` and `.im` components.
    @variables x::Complex
    @test !isnothing(metadata(unwrap(x)))
    @test_nowarn substitute(x + 1.0, x => 1.0im)
end

@testset "getname" begin
    @variables t a b x::Complex y(t)::Complex z(a, b)::Complex
    @test hasname(x)
    @test getname(x) == :x
    @test hasname(y)
    @test getname(y) == :y
    @test hasname(z)
    @test getname(z) == :z
    @test !hasname(2x)
    @test !hasname(x + y)
end
