using Test
using Symbolics

@testset "historical complex PR regressions" begin
    @testset "#1492 linear expansion stays generic" begin
        @variables x::Real y::Real

        a, b, islinear = Symbolics.linear_expansion(2im * x + im, x)
        @test islinear
        @test isequal(a, 2im)
        @test isequal(b, im)

        a, b, islinear = Symbolics.linear_expansion(im * x + im * y, x)
        @test islinear
        @test isequal(a, im)
        @test isequal(b, im * y)

        _, _, islinear = Symbolics.linear_expansion(im * x^2 + im * x, x)
        @test !islinear
    end

    @testset "#160 tuple-argument build_function" begin
        @variables a::Real b::Real
        out = a + im * b
        f = build_function(out, (a, b); expression = Val(false))
        @test f((1.0, 2.0)) == 1.0 + 2.0im
    end
end
