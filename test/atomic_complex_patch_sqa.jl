using Pkg

"""
Patch a temporary installed/developed SecondQuantizedAlgebra checkout just enough to let
its current code run against the atomic Symbolics prototype. This is test scaffolding, not
an upstream SQA change: it deliberately removes only representation assumptions that the
new Symbolics API is meant to obsolete.
"""
function patch_sqa_for_atomic_symbolics!()
    info = only(v for v in values(Pkg.dependencies()) if v.name == "SecondQuantizedAlgebra")
    src = info.source
    src === nothing && error("SecondQuantizedAlgebra source path unavailable")

    # SQA currently uses Symbolics.IM solely to avoid the eager Complex{Num} path. With an
    # atomic numeric wrapper, Base.im is the intended literal coefficient.
    for (root, _, files) in walkdir(joinpath(src, "src"))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(root, file)
            text = read(path, String)
            occursin("Symbolics.IM", text) || continue
            write(path, replace(text, "Symbolics.IM" => "im"))
        end
    end

    # Current SQA's prefactor printer only recognizes Complex{Num}. Let the temporary
    # checkout display an atomic SymbolicNumber directly so precompilation can proceed to
    # the algebraic tests. The eventual SQA migration should generalize this API properly.
    printing = joinpath(src, "src", "printing", "printing.jl")
    open(printing, "a") do io
        write(io, "\nif isdefined(Symbolics, :SymbolicNumber)\n")
        write(io, "    show_display(io::IO, c::Symbolics.SymbolicNumber) = print(io, c)\n")
        write(io, "end\n")
    end

    return src
end
