using Test
using Symbolics
using SymbolicUtils

@testset "atomic complex arithmetic without constructor interception" begin
    @variables x::Real y::Real z::Complex

    # Explicit Cartesian construction keeps the Base representation contract.
    cart = Complex(x, y)
    @test cart isa Complex{Num}
    @test isequal(real(cart), x)
    @test isequal(imag(cart), y)

    # Ordinary complex arithmetic never routes through that Cartesian constructor.
    for ex in (
        im * x,
        x * im,
        x + 3im,
        3im + x,
        x - 3im,
        3im - x,
        x * (2 + 3im),
        (2 + 3im) * x,
        x / (2 + 3im),
        (2 + 3im) / x,
    )
        @test ex isa Symbolics.SymbolicNumber
        @test !(ex isa Complex{Num})
        @test SymbolicUtils.symtype(Symbolics.unwrap(ex)) <: Number
    end

    # A declared complex scalar remains atomic independently of literal-complex mixing.
    @test z isa Symbolics.SymbolicNumber
    @test exp(z) isa Symbolics.SymbolicNumber
    @test z + x isa Symbolics.SymbolicNumber

    f = Symbolics.build_function(x + 3im, x; expression = Val(false))
    @test f(2.0) == 2 + 3im
end
