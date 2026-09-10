using Test
using Symbolics
using SymbolicUtils

@testset "numeric wrapper lattice" begin
    r = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_real; type = Real)
    n = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_number; type = Number)
    c = SymbolicUtils.Sym{Symbolics.VartypeT}(:wrapper_complex; type = Complex{Real})

    @test Symbolics.wrapper_type(Real) === Num
    @test Symbolics.wrapper_type(Number) === Symbolics.SymbolicNumber
    @test Symbolics.wrapper_type(Complex{Real}) === Symbolics.SymbolicNumber

    @test Symbolics.wrap(r) isa Num
    @test Symbolics.wrap(n) isa Symbolics.SymbolicNumber
    @test Symbolics.wrap(c) isa Symbolics.SymbolicNumber

    # The wide wrapper may contain a narrower numeric symtype when container/promotion
    # stability requires it. Re-running `wrap` on the raw node recovers the narrow wrapper.
    wide_zero = zero(Symbolics.SymbolicNumber)
    @test wide_zero isa Symbolics.SymbolicNumber
    @test SymbolicUtils.symtype(Symbolics.unwrap(wide_zero)) <: Real
    @test Symbolics.wrap(Symbolics.unwrap(wide_zero)) isa Num

    mixed = [Symbolics.wrap(r), Symbolics.wrap(c)]
    @test eltype(mixed) === Symbolics.SymbolicNumber
    @test mixed[1] isa Symbolics.SymbolicNumber
    @test SymbolicUtils.symtype(Symbolics.unwrap(mixed[1])) <: Real
    @test SymbolicUtils.symtype(Symbolics.unwrap(mixed[2])) <: Complex
end

# A package-defined numeric symbolic domain must remain able to register a more-specific
# wrapper. The built-in Number fallback must not steal it.
abstract type ProbeNumericDomain <: Number end
@symbolic_wrap struct ProbeNumericWrapper <: ProbeNumericDomain
    val::SymbolicUtils.BasicSymbolic{Symbolics.VartypeT}
end
SymbolicUtils.unwrap(x::ProbeNumericWrapper) = x.val

@testset "custom numeric wrapper specificity" begin
    p = SymbolicUtils.Sym{Symbolics.VartypeT}(:probe_numeric; type = ProbeNumericDomain)
    @test Symbolics.wrapper_type(ProbeNumericDomain) === ProbeNumericWrapper
    @test Symbolics.wrap(p) isa ProbeNumericWrapper
    @test Symbolics.unwrap(Symbolics.wrap(p)) === p

    # Existing built-ins retain their more-specific registrations.
    @test Symbolics.wrapper_type(Real) === Num
    @test Symbolics.wrapper_type(Complex{Real}) === Symbolics.SymbolicNumber
end
