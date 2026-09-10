using Pkg

"""
Patch a temporary SecondQuantizedAlgebra checkout just enough to test it against atomic
Symbolics. These are downstream migration shims, not upstream Symbolics behavior.
"""
function patch_sqa_for_atomic_symbolics!(src::AbstractString)
    # `Symbolics.IM` existed only to avoid the old eager `Complex{Num}` path. Atomic
    # symbolic numbers make the literal Base.im coefficient the correct representation.
    for (root, _, files) in walkdir(joinpath(src, "src"))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(root, file)
            text = read(path, String)
            occursin("Symbolics.IM", text) || continue
            write(path, replace(text, "Symbolics.IM" => "im"))
        end
    end

    # SQA's current prefactor printer assumes every non-native symbolic coefficient lowers
    # to `Complex{Num}`. Under atomic Symbolics it may lower to `Num` or SymbolicNumber.
    printing = joinpath(src, "src", "printing", "printing.jl")
    open(printing, "a") do io
        write(io, "\nshow_display(io::IO, c::Num) = print(io, c)\n")
        write(io, "needs_pf_parens(c::Num) = is_loose_head(c)\n")
        write(io, "if isdefined(Symbolics, :SymbolicNumber)\n")
        write(io, "    show_display(io::IO, c::Symbolics.SymbolicNumber) = print(io, c)\n")
        write(io, "    needs_pf_parens(c::Symbolics.SymbolicNumber) = iszero(imag(c)) && is_loose_head(real(c))\n")
        write(io, "end\n")
    end

    return src
end

function patch_sqa_for_atomic_symbolics!()
    info = only(v for v in values(Pkg.dependencies()) if v.name == "SecondQuantizedAlgebra")
    src = info.source
    src === nothing && error("SecondQuantizedAlgebra source path unavailable")
    return patch_sqa_for_atomic_symbolics!(src)
end
