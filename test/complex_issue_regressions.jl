using Test
using LinearAlgebra
using SpecialFunctions
using Symbolics
using SymbolicUtils

const SN = Symbolics.SymbolicNumber
raw(x) = Symbolics.unwrap(x)

@testset "historical complex issue regressions" begin
    @testset "#118 #1674 #718 #883 #1830 elementary complex functions" begin
        @variables x::Real z::Complex

        for ex in (exp(im * x), sqrt(-im + x), log(-im + x), sin(z), cos(z), exp(z), sqrt(z), log(z))
            @test ex isa SN
            @test !(ex isa Complex{Num})
        end
        @test SymbolicUtils.operation(raw(exp(im * x))) === exp
    end

    @testset "#232 derivative of exp(im*x)" begin
        @variables x::Real
        dex = expand_derivatives(Differential(x)(exp(im * x)))
        f = build_function(dex, x; expression = Val(false))
        @test f(0.37) ≈ im * exp(0.37im)
    end

    @testset "#327 #921 variable discovery" begin
        @variables t::Real x::Real u::Real v::Real z::Complex
        explicit = Complex(u, v)
        vars_atomic = Set(Symbolics.get_variables(x + t * z + x))
        vars_cart = Set(Symbolics.get_variables(x + t * explicit + x))
        @test raw(t) in vars_atomic
        @test raw(x) in vars_atomic
        @test raw(z) in vars_atomic
        @test raw(t) in vars_cart
        @test raw(x) in vars_cart
        @test raw(u) in vars_cart
        @test raw(v) in vars_cart
    end

    @testset "#534 #905 #1109 #1813 substitution and domains" begin
        @variables x::Real z::Number f::Real a::Real b::Real c::Real

        @test Symbolics.value(substitute(im * z, Dict(z => im); fold = Val(true))) == -1
        p = 0.4 + 1.7im * z
        @test Symbolics.value(substitute(p, Dict(z => 0.2 + 1.0im); fold = Val(true))) ≈ 0.4 + 1.7im * (0.2 + 1.0im)

        L = a * z^2 + b * z + c
        subL = substitute(L, Dict(z => 2pi * f * im); fold = Val(false))
        @test subL isa Union{SN, Num, SymbolicUtils.BasicSymbolic}
        Lf = build_function(subL, a, b, c, f; expression = Val(false))
        @test Lf(2.0, 3.0, 4.0, 0.7) ≈ 2.0 * (2pi * 0.7im)^2 + 3.0 * (2pi * 0.7im) + 4.0

        # A variable declared Real must not silently become a semantically-complex Num.
        @test_throws Exception substitute(x + 1, Dict(x => 1 + 2im); fold = Val(false))
    end

    @testset "#159 #311 #354 #1016 #1391 build_function and expand" begin
        @variables a::Real b::Real u::Complex k1::Real k2::Real

        explicit = Complex(a, b)
        ef = build_function(explicit, a, b; expression = Val(false))
        @test ef(2.0, -3.0) == 2.0 - 3.0im

        uf = build_function(u^2 + exp(u), u; expression = Val(false))
        @test uf(1.0 + 0.5im) ≈ (1.0 + 0.5im)^2 + exp(1.0 + 0.5im)

        A = [u, im * u, conj(u)]
        af = build_function(A, u; expression = Val(false))[1]
        @test af(1 + 2im) == [1 + 2im, -2 + 1im, 1 - 2im]

        elem = exp(pi * im / 4) * cos(k1) + exp(-pi * im / 4) * cos(k2)
        expanded = expand(elem)
        @test expanded isa SN
        bf = build_function(expanded, k1, k2; expression = Val(false))
        @test bf(0.2, 0.4) ≈ exp(pi * im / 4) * cos(0.2) + exp(-pi * im / 4) * cos(0.4)

        @variables q::Real
        df = expand_derivatives(Differential(q)(1.0 + exp(im * q)))
        dfun = build_function(df, q; expression = Val(false))
        @test dfun(0.0) ≈ im
    end

    @testset "#341 expansion and #1116 degree" begin
        @variables x::Real y::Real z::Real
        ex = (x + im * y)^4
        expanded = expand(ex)
        f = build_function(expanded, x, y; expression = Val(false))
        @test f(0.3, -0.7) ≈ (0.3 - 0.7im)^4
        @test Symbolics.degree(im + z, z) == 1
    end

    @testset "#777 rational construction semantics" begin
        @variables x::Real
        # `//` constructs exact rationals and is intentionally restricted to rational
        # domains. A generic symbolic Real/complex expression must use algebraic division.
        for term in (im + x, 1 + im + x)
            @test_throws MethodError term // term
            r = term / term
            rf = build_function(r, x; expression = Val(false))
            @test rf(0.4) ≈ 1
        end
    end

    @testset "#800 printing" begin
        @variables z::Complex
        s = sprint(show, z + im)
        @test occursin("z", s)
        @test occursin("im", s)
        @test !occursin("real(", s)
        @test !occursin("imag(", s)
    end

    @testset "#832 real/imag simplification" begin
        @variables r1::Real r2::Real i1::Real i2::Real
        x1 = r1 + i1 * im
        x2 = r2 + i2 * im
        got = simplify(real(x1 * x2); expand = true)
        expected = r1 * r2 - i1 * i2
        gf = build_function(got - expected, r1, r2, i1, i2; expression = Val(false))
        @test gf(1.2, -0.4, 0.7, 2.0) ≈ 0
    end

    @testset "#861 complex symbolic LinearAlgebra" begin
        @variables Ω::Real ω0::Real Δ::Real
        M = [-ω0 2im * Ω 0; -2im * Ω -ω0 2im * Δ; 0 -2im * Δ -ω0]
        d = det(M)
        df = build_function(d, Ω, ω0, Δ; expression = Val(false))
        vals = (0.7, 1.3, -0.2)
        Mnum = [-vals[2] 2im * vals[1] 0; -2im * vals[1] -vals[2] 2im * vals[3]; 0 -2im * vals[3] -vals[2]]
        @test df(vals...) ≈ det(Mnum)
    end

    @testset "#884 compact atomic rational expression" begin
        @variables z::Complex
        ex = 1 / (1 - z^10)
        @test ex isa SN
        @test length(sprint(show, ex)) < 200
        f = build_function(ex, z; expression = Val(false))
        @test f(0.3 + 0.2im) ≈ 1 / (1 - (0.3 + 0.2im)^10)
    end

    @testset "#894 sound complex differential expression" begin
        @variables x::Real z(x)::Complex
        D = Differential(x)
        ex = x * D(z) + z
        @test ex isa SN
        @test !(ex isa Num)
        eq = D(z) ~ ex
        @test eq isa Equation
    end

    @testset "#1199 #1485 complex equations and linear solve" begin
        @variables x::Real
        eq = x + 3 + im ~ 0
        @test eq isa Equation
        f = build_function(eq.lhs, x; expression = Val(false))
        @test f(2.0) == 5 + im
        @test solve_for(x + im, x) == -im
    end

    @testset "#1487 left division" begin
        @variables x::Complex y::Complex
        ex = x \ y
        f = build_function(ex, x, y; expression = Val(false))
        @test f(1 + 2im, 3 - im) ≈ ((1 + 2im) \ (3 - im))
    end

    @testset "#1661 heterogeneous numeric symbolic domains" begin
        @variables t::Real x::Number y::Complex z(t)::Real
        v = [t, x, y, z]
        @test eltype(v) == SN
        @test length(v) == 4
    end

    @testset "#1372 dot follows Julia Hermitian semantics" begin
        @variables mass::Real qsqu::Real
        q2 = [0, 0, -(mass^2 + qsqu) / sqrt(qsqu) / 2,
              -im * (mass^2 + qsqu) / sqrt(qsqu) / 2]
        got = simplify(dot(q2, q2))
        direct = simplify(sum(q2[i] * q2[i] for i in eachindex(q2)))
        gdot = build_function(got, mass, qsqu; expression = Val(false))
        gdirect = build_function(direct, mass, qsqu; expression = Val(false))
        m, q = 1.4, 2.3
        qnum = [0, 0, -(m^2 + q) / sqrt(q) / 2, -im * (m^2 + q) / sqrt(q) / 2]
        @test gdot(m, q) ≈ dot(qnum, qnum)
        @test gdirect(m, q) ≈ sum(v * v for v in qnum)
        @test !isapprox(gdot(m, q), gdirect(m, q))
    end

    @testset "#465 complex symbolic arrays scalarize" begin
        @variables (v::Complex)[1:2]
        got = scalarize([1 2; 3 4] * v)
        @test length(got) == 2
        @test isequal(got[1], v[1] + 2v[2])
        @test isequal(got[2], 3v[1] + 4v[2])
    end

    @testset "#645 symbolic arrays and complex substitution" begin
        @variables r[1:1, 1:1]::Real phi[1:1, 1:1]::Real
        H = hankelh1.(0, r)
        Q = im .* cos.(phi)
        y = scalarize((Q .* H)[1])
        sub = substitute(y, Dict(r[1, 1] => 2.0, phi[1, 1] => 0.3); fold = Val(true))
        @test Symbolics.value(sub) ≈ im * cos(0.3) * hankelh1(0, 2.0)
    end

    @testset "#577 do not implicitly split complex equations" begin
        @variables z::Complex w::Complex
        eq = z ~ w
        @test eq isa Equation
        @test !(eq isa AbstractArray)
        f = build_function(eq.lhs - eq.rhs, z, w; expression = Val(false))
        @test f(1 + 2im, 0.5 - im) ≈ 0.5 + 3im
    end

    @testset "#558 #1011 complex differentiation contract" begin
        @variables t::Real w(t)::Complex
        D = Differential(t)
        dw = expand_derivatives(D(w))
        @test dw != 0
        @test isequal(expand_derivatives(D(conj(w))), conj(dw))
        @test isequal(expand_derivatives(D(real(w))), real(dw))
        @test isequal(expand_derivatives(D(imag(w))), imag(dw))

        @variables z::Complex
        Dz = Differential(z)
        @test expand_derivatives(Dz(z)) == 1
        # Non-holomorphic projections remain unevaluated rather than silently claiming zero.
        @test Symbolics.is_derivative(expand_derivatives(Dz(conj(z))))
        @test Symbolics.is_derivative(expand_derivatives(Dz(real(z))))
        @test Symbolics.is_derivative(expand_derivatives(Dz(imag(z))))
    end
end
