using Test
using Symbolics
using SymbolicUtils
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using SciMLBase

@testset "ModelingToolkit atomic numeric state" begin
    @variables x(t)::Number
    @test x isa Symbolics.SymbolicNumber

    @mtkcompile sys = System([D(x) ~ -conj(x)], t)
    prob = ODEProblem(sys, [x => 1.0 + 3im], (0.0, 0.2))
    sol = solve(prob, Tsit5())
    @test SciMLBase.successful_retcode(sol)
    @test sol.u[1][1] == 1.0 + 3im
    @test eltype(sol.u[1]) <: Complex
end

@testset "ModelingToolkit explicit complex domain admission" begin
    # Current ModelingToolkitBase accepts `Complex` exactly but accidentally rejects
    # concrete complex domains such as `Complex{Real}`. Record that downstream gap, then
    # apply the minimal proposed fix so the rest of the end-to-end semantics are tested.
    @test !ModelingToolkitBase.is_floatingpoint_symtype(Complex{Real})
    @eval ModelingToolkitBase is_floatingpoint_symtype(::Type{T}) where {T <: Complex} = true
    @test ModelingToolkitBase.is_floatingpoint_symtype(Complex{Real})

    @variables z(t)::Complex{Real}
    @test z isa Symbolics.SymbolicNumber
    @test SymbolicUtils.symtype(Symbolics.unwrap(z)) <: Complex

    rhs = im * z
    @test rhs isa Symbolics.SymbolicNumber
    @test !(rhs isa Complex{Num})

    @mtkcompile sysc = System([D(z) ~ im * z], t)
    probc = ODEProblem(sysc, [z => 1.0 + 0im], (0.0, 0.2))
    solc = solve(probc, Tsit5(), saveat = 0.2)
    @test SciMLBase.successful_retcode(solc)
    @test solc.u[end][1] ≈ exp(0.2im) rtol = 1e-5
end
