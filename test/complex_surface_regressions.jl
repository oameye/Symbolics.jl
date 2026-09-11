using Test
using Symbolics
using SymbolicUtils
using Latexify
using LinearAlgebra
using SparseArrays

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

    @testset "general numeric wrapper reaches high-level differentiation" begin
        @variables z::Complex w::Complex

        dz = Symbolics.derivative(z^2 + im * z, z)
        @test dz isa SN
        @test iszero(simplify(dz - (2z + im)))

        g = Symbolics.gradient(z * w + im * z, [z, w])
        @test length(g) == 2
        @test iszero(simplify(g[1] - (w + im)))
        @test iszero(simplify(g[2] - z))

        J = Symbolics.jacobian([z^2 + w, im * z + w^2], [z, w])
        @test size(J) == (2, 2)
        @test iszero(simplify(J[1, 1] - 2z))
        @test isone(simplify(J[1, 2]))
        @test iszero(simplify(J[2, 1] - im))
        @test iszero(simplify(J[2, 2] - 2w))

        Js = Symbolics.sparsejacobian([z^2 + w, im * z + w^2], [z, w])
        @test Js isa SparseMatrixCSC
        @test iszero(simplify(Js[1, 1] - 2z))
        @test iszero(simplify(Js[2, 1] - im))

        H = Symbolics.hessian(im * z^2 + z * w, [z, w])
        @test size(H) == (2, 2)
        @test iszero(simplify(H[1, 1] - 2im))
        @test isone(simplify(H[1, 2]))
        @test isone(simplify(H[2, 1]))
        @test iszero(simplify(H[2, 2]))

        Hs = Symbolics.sparsehessian(im * z^2 + z * w, [z, w])
        @test Hs isa SparseMatrixCSC
        @test iszero(simplify(Hs[1, 1] - 2im))
        @test isone(simplify(Hs[1, 2]))
        @test isone(simplify(Hs[2, 1]))
    end

    @testset "general numeric wrapper reaches symbolic linear algebra" begin
        @variables z::Complex w::Complex
        M = [z 1; 1 w]

        @test lu(M; check = false) isa LinearAlgebra.LU
        @test iszero(simplify(det(M; laplace = false) - det(M; laplace = true)))

        Minv = inv(M; laplace = false)
        ident = simplify.(M * Minv)
        @test isone(ident[1, 1])
        @test iszero(ident[1, 2])
        @test iszero(ident[2, 1])
        @test isone(ident[2, 2])

        Mex = exp([z zero(z); zero(w) w])
        @test Mex isa Symbolics.Arr{SN, 2}
    end

    @testset "complex symbolic arrays preserve result domains" begin
        @variables (v::Complex)[1:2]
        nv = norm(v)
        @test nv isa Num
        @test SymbolicUtils.symtype(Symbolics.unwrap(nv)) <: Real
    end

    @testset "complex symbolic linear systems" begin
        @variables z::Complex w::Complex

        scalar_sol = symbolic_linear_solve(z + im ~ 0, z)
        @test iszero(simplify(scalar_sol + im))

        sols = symbolic_linear_solve([z + w ~ 1, z - w ~ im], [z, w])
        @test length(sols) == 2
        @test iszero(simplify(sols[1] - (1 + im) / 2))
        @test iszero(simplify(sols[2] - (1 - im) / 2))
    end

    @testset "semi-polynomial forms preserve complex coefficients" begin
        @variables x::Real y::Real
        expr = im * x + (1 + im) * y + 2

        A, c = semilinear_form([expr], [x, y])
        @test iszero(simplify(A[1, 1] - im))
        @test iszero(simplify(A[1, 2] - (1 + im)))
        @test iszero(simplify((A * [x, y] + c)[1] - expr))

        qexpr = im * x^2 + (1 + im) * x * y + 2y + 3
        Aq, Bq, v2, cq = semiquadratic_form([qexpr], [x, y])
        @test iszero(simplify((Aq * [x, y] + Bq * v2 + cq)[1] - qexpr))
    end
end
