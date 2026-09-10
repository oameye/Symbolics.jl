using Test
using Symbolics
using SymbolicUtils
using Latexify

const SN = Symbolics.SymbolicNumber

@testset "complex feature surfaces" begin
    @testset "#1391 cis construction, differentiation, and codegen" begin
        @variables x::Real

        c = cis(x)
        @test c isa SN
        @test !(c isa Complex{Num})

        fexpr = 1 + c
        f = build_function(fexpr, x; expression = Val(false))
        @test f(0.37) ≈ 1 + cis(0.37)

        D = Differential(x)
        dex = expand_derivatives(D(fexpr))
        df = build_function(dex, x; expression = Val(false))
        @test df(0.37) ≈ im * cis(0.37)
    end

    @testset "#1674 compact conjugated complex exponential differentiation" begin
        @variables y::Real
        D = Differential(y)
        phase = exp(im * y)
        @test phase isa SN
        dex = expand_derivatives(D(conj(phase)))
        @test !Symbolics.is_derivative(dex)
        f = build_function(dex, y; expression = Val(false))
        @test f(0.37) ≈ conj(im * exp(0.37im))
    end

    @testset "#389 #416 complex Latexify" begin
        # Exact raw SymbolicUtils-style MWE from #416.
        @syms sx::Real
        raw_tex = latexify(im * sx)
        @test !isempty(string(raw_tex))

        @variables x::Real z::Complex a::Real b::Real
        expressions = (z, im * x, exp(im * x), Complex(a, b))
        for ex in expressions
            tex = latexify(ex)
            @test !isempty(string(tex))
        end

        ztex = sprint(show, MIME"text/latex"(), z)
        @test !isempty(ztex)
        @test !occursin("real(z)", ztex)
        @test !occursin("imag(z)", ztex)

        phasetex = sprint(show, MIME"text/latex"(), exp(im * x))
        @test !isempty(phasetex)
        @test occursin("exp", phasetex) || occursin("e", lowercase(phasetex))
    end
end
