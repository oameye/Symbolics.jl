using Test
using Symbolics
using SymbolicUtils
using SecondQuantizedAlgebra
using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, mtkcompile, System

@testset "QuantumCumulants against atomic Symbolics" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables Δ::Real η::Real κ::Real

    H = Δ * a' * a + η * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ])

    mmap = moment_variable_map(eqs)
    u_raw = first(values(mmap))
    @test SymbolicUtils.symtype(SymbolicUtils.unwrap(u_raw)) == Number
    @test Symbolics.wrap(SymbolicUtils.unwrap(u_raw)) isa Symbolics.SymbolicNumber

    sys = System(eqs; name = :atomic_complex_qc)
    sysc = mtkcompile(sys)
    @test !isempty(unknowns(sysc))

    uw = Symbolics.wrap(SymbolicUtils.unwrap(first(unknowns(sysc))))
    if SymbolicUtils.symtype(Symbolics.unwrap(uw)) <: Number
        @test !isequal(conj(uw), uw) || SymbolicUtils.symtype(Symbolics.unwrap(uw)) <: Real
    end
end
