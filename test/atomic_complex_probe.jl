using Test
using Symbolics
using SymbolicUtils

@testset "atomic complex scalar probe" begin
    @variables x::Real z::Complex n::Number

    @test x isa Num
    @test z isa Symbolics.SymbolicNumber
    @test n isa Symbolics.SymbolicNumber
    @test SymbolicUtils.symtype(Symbolics.unwrap(z)) <: Complex
    @test SymbolicUtils.symtype(Symbolics.unwrap(n)) == Number

    zx = im * x
    @test zx isa Symbolics.SymbolicNumber
    @test Set(Symbolics.get_variables(zx)) == Set([Symbolics.unwrap(x)])
    @test !isdefined(Symbolics, :IM)

    @testset "complex projections and conjugation" begin
        rz = real(z)
        iz = imag(z)
        cz = conj(z)
        @test rz isa Num
        @test iz isa Num
        @test cz isa Symbolics.SymbolicNumber
        @test !isequal(cz, z)
        @test SymbolicUtils.operation(Symbolics.unwrap(cz)) === conj
    end

    @testset "elementary functions remain atomic" begin
        for f in (exp, sin, cos, log, sqrt)
            y = f(z)
            @test y isa Symbolics.SymbolicNumber
            @test !(y isa Complex{Num})
            @test SymbolicUtils.operation(Symbolics.unwrap(y)) === f
        end

        phase = exp(im * x)
        @test phase isa Symbolics.SymbolicNumber
        @test !(phase isa Complex{Num})
        @test SymbolicUtils.operation(Symbolics.unwrap(phase)) === exp

        # Regression shape for historical sqrt/log/exp Complex{Num} failures: the
        # non-real argument remains one symbolic expression and never enters Base's
        # numerical Complex algorithms during construction.
        for f in (sqrt, log, exp)
            y = f(-im + x)
            @test y isa Symbolics.SymbolicNumber
            @test !(y isa Complex{Num})
            @test SymbolicUtils.operation(Symbolics.unwrap(y)) === f
        end
    end

    @testset "substitution and code generation" begin
        @test Symbolics.value(substitute(im * x, Dict(x => 2.0); fold = Val(true))) == 2.0im
        @test Symbolics.value(substitute(z, Dict(z => 1.0 + 2.0im); fold = Val(true))) == 1.0 + 2.0im
        @test Symbolics.value(substitute(im * z, Dict(z => 1.0 + 2.0im); fold = Val(true))) == -2.0 + 1.0im

        f = Symbolics.build_function(im * x, x; expression = Val(false))
        @test f(2.0) == 2.0im

        g = Symbolics.build_function(exp(im * x), x; expression = Val(false))
        @test g(0.25) ≈ exp(0.25im)
    end

    @testset "real-variable differentiation representation" begin
        @variables t::Real w(t)::Complex
        D = Differential(t)
        dw = D(w)
        @test dw isa Symbolics.SymbolicNumber
        @test SymbolicUtils.symtype(Symbolics.unwrap(dw)) <: Number
        @test SymbolicUtils.operation(Symbolics.unwrap(dw)) isa Differential

        # Do not impose a holomorphic/non-holomorphic calculus rule here. The atomic
        # representation must preserve both expressions without splitting into re/im;
        # their mathematical relationship is tested in the dedicated calculus probe.
        dcw = D(conj(w))
        @test dcw isa Symbolics.SymbolicNumber
        @test SymbolicUtils.operation(Symbolics.unwrap(dcw)) isa Differential
    end

    @testset "promotion and small arrays" begin
        @variables z1::Complex z2::Complex
        a = z1 + 2z2
        b = 3.0 * z1 - im * z2
        @test a isa Symbolics.SymbolicNumber
        @test b isa Symbolics.SymbolicNumber

        A = [z1 z2; conj(z1) z1 + z2]
        @test eltype(A) <: Number
        h, h! = Symbolics.build_function(A, z1, z2; expression = Val(false))
        out = h(1.0 + 2.0im, 3.0 - 1.0im)
        @test out == [1.0 + 2.0im 3.0 - 1.0im; 1.0 - 2.0im 4.0 + 1.0im]

        dest = similar(out)
        h!(dest, 1.0 + 2.0im, 3.0 - 1.0im)
        @test dest == out
    end

    @testset "historical representation regressions" begin
        # Heterogeneous numeric symbolic domains must have a common atomic wrapper rather
        # than forcing every element through Complex{Num} conversion.
        @variables t0::Real generic::Number cplx::Complex q(t0)::Real
        v = [t0, generic, cplx, q]
        @test eltype(v) == Symbolics.SymbolicNumber
        @test length(v) == 4

        # Equation deliberately stores raw BasicSymbolic payloads. Test semantic
        # preservation through code generation instead of requiring a wrapper there.
        eq = x + 3 + im ~ 0
        @test eq.lhs isa SymbolicUtils.BasicSymbolic
        eqf = Symbolics.build_function(eq.lhs, x; expression = Val(false))
        @test eqf(0.0) == 3 + im

        rational_complex = 1 / (1 - z^10)
        @test rational_complex isa Symbolics.SymbolicNumber
        @test !(rational_complex isa Complex{Num})
        @test SymbolicUtils.operation(Symbolics.unwrap(z^10)) === ^
    end
end
