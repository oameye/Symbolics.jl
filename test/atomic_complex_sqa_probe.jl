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

    # Upstream should make a literal imaginary coefficient safe: no fake Symbolics.IM
    # should be required merely to keep the average complex-valued.
    avg = average((2 + 3im) * a)
    @test SymbolicUtils.symtype(avg) <: Number

    # The ordinary upstream exponential now has the compact representation SQA's expim
    # was originally introduced to protect from Complex{Num} expansion.
    upstream_phase = exp(im * ω * t)
    @test upstream_phase isa Symbolics.SymbolicNumber
    @test !(upstream_phase isa Complex{Num})

    sqa_phase = expim(ω * t)
    @test !iszero(sqa_phase)
end
