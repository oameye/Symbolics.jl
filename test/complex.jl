using Symbolics, Test
using SymbolicUtils: metadata, symtype, unwrap_const
using Symbolics: unwrap
using SymbolicIndexingInterface: getname, hasname

@variables a b::Real z::Complex (Z::Complex)[1:10]

@testset "legacy Complex{Num} contracts" begin
    @test a isa Num
    @test b isa Num
    @test eltype(Z) <: Complex{Num}

    for x in [z, Z[1], z+a, z*a, z^2, z/z] # z/z is sus
        @test x isa Complex{Num}
        @test real(x) isa Num
        @test imag(x) isa Num
        @test conj(x) isa Complex{Num}
    end

    # issue #314
    bi = a+a*im
    bs = substitute(bi, (Dict(a=>1.0))) # returns 1.0 + im
    @test bs isa Complex{Num}
    bv = unwrap_const(Symbolics.value(bs))
    @test typeof(bv) == ComplexF64
end

@testset "legacy repr" begin
    @test repr(z) == "z"
    @test repr(a + b*im) == "a + b*im"
end

@testset "legacy metadata" begin
    z1 = z+1.0
    @test_nowarn substitute(z1, z=>1.0im)
    @test metadata(z1) == unwrap(z1.im).metadata
    @test metadata(z1) == unwrap(z1.re).metadata
    z2 = 1.0 + z*im
    @test isnothing(metadata(unwrap(z1.re)))
end

@testset "legacy getname" begin
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

@testset "atomic complex scalar invariants" begin
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

    # `im` is a coefficient, never a free symbolic variable.
    @test Symbolics.get_variables(ai) == [a]
    @test Symbolics.get_variables(ia) == [a]

    @test repr(ai) == "a*im" || repr(ai) == "im*a"
    @test repr(ia) == "a*im" || repr(ia) == "im*a"

    @test Symbolics.value(substitute(ai, Dict(a => 2.0))) == 2.0im
    @test Symbolics.value(substitute(ia, Dict(a => 2.0))) == 2.0im

    f = Symbolics.build_function(ai, a; expression = Val(false))
    @test f(2.0) == 2.0im
end

@testset "atomic substitution" begin
    bi = a + a * im
    @test bi isa Symbolics.SymbolicNumber

    bs = substitute(bi, Dict(a => 1.0))
    @test bs isa Symbolics.SymbolicNumber
    bv = unwrap_const(Symbolics.value(bs))
    @test bv == 1.0 + 1.0im
    @test typeof(bv) == ComplexF64
end

@testset "compact atomic representation" begin
    @test repr(z) == "z"
    @test repr(a + b * im) == "a + b*im"
    @test repr(exp(im * a)) == "exp(im*a)" || repr(exp(im * a)) == "exp(a*im)"
end

@testset "atomic metadata" begin
    @variables x::Complex
    @test !isnothing(metadata(unwrap(x)))
    @test_nowarn substitute(x + 1.0, x => 1.0im)
end

@testset "atomic getname" begin
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
