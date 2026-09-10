using Test
using Symbolics
using SymbolicUtils
using SecondQuantizedAlgebra
import SecondQuantizedAlgebra: expim

@testset "SecondQuantizedAlgebra against atomic Symbolics" begin
    h = FockSpace(:f)
    @qnumbers a::Destroy(h)
    @variables t::Real ω::Real

    leaf = average(a)
    @test SymbolicUtils.symtype(leaf) == Number

    lifted = make_time_dependent(leaf, t)
    @test SymbolicUtils.symtype(lifted) == Number
    @test Symbolics.wrap(lifted) isa Symbolics.SymbolicNumber

    avg = average((2 + 3im) * a)
    @test SymbolicUtils.symtype(avg) <: Number

    upstream_phase = exp(im * ω * t)
    @test upstream_phase isa Symbolics.SymbolicNumber
    @test !(upstream_phase isa Complex{Num})

    sqa_phase = expim(ω * t)
    @test !iszero(sqa_phase)

    # Exercise the prefactor display path that previously assumed `Complex{Num}`.
    @test !isempty(sprint(show, avg))
end
