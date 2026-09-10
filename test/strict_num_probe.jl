using Test
using Symbolics
using SymbolicUtils
using LinearAlgebra

@testset "strict Num domain" begin
    @variables x::Real y::Real z::Complex

    # Operations closed over the reals remain represented by Num.
    for ex in (x + y, x * y, sin(x), cos(x), exp(x), abs(x), x^2)
        @test ex isa Num
        @test SymbolicUtils.symtype(Symbolics.unwrap(ex)) <: Real
    end

    # Operations whose mathematical codomain may leave R must never be hidden in Num.
    for ex in (sqrt(x), log(x), x^(1//2))
        @test !(ex isa Num) || SymbolicUtils.symtype(Symbolics.unwrap(ex)) <: Real
    end

    # Literal imaginary arithmetic takes the atomic non-real wrapper path.
    @test im * x isa Symbolics.SymbolicNumber
    @test x + im * y isa Symbolics.SymbolicNumber
    @test z + x isa Symbolics.SymbolicNumber

    @testset "real linear algebra" begin
        A = [x 1; y 2]
        B = [2 3; 4 5]
        C = A * B
        @test size(C) == (2, 2)
        @test all(c -> c isa Num, C)

        # Exercise symbolic LU storage and ordinary-number insertion paths mentioned
        # in Num's historical relaxed-domain comment.
        M = [x + 1 one(x); one(x) y + 2]
        F = lu(M; check = false)
        @test size(F.L) == (2, 2)
        @test size(F.U) == (2, 2)
    end

    @testset "mixed/complex linear algebra" begin
        A = [z one(z); x im * y]
        @test eltype(A) <: Number
        C = A * [1 2; 3 4]
        @test size(C) == (2, 2)
        @test all(c -> c isa Symbolics.SymbolicNumber, C)

        f_oop, f_iip = Symbolics.build_function(C, z, x, y; expression = Val(false))
        expected = [4 + 2im 6 + 4im; 2 + 3im 4 + 4im]
        @test f_oop(1 + 2im, 2.0, 1.0) == expected
        out = similar(expected)
        f_iip(out, 1 + 2im, 2.0, 1.0)
        @test out == expected
    end
end
