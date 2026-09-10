"""
Patch a temporary QuantumCumulants checkout to probe the atomic symbolic-scalar boundary.
The final QC migration should retain a finite concrete wrapper union; this test shim uses
`Number` only to discover every place that still assumes `Num` can contain complex trees.
"""
function patch_qc_for_atomic_symbolics!(src::AbstractString)
    for (root, _, files) in walkdir(joinpath(src, "src"))
        for file in files
            endswith(file, ".jl") || continue
            path = joinpath(root, file)
            text = read(path, String)
            text = replace(text, "Symbolics.Num(" => "Symbolics.wrap(")
            if endswith(path, joinpath("src", "moments.jl"))
                text = replace(text,
                    "_reduce_ground_in_drift(x::Symbolics.Num)" =>
                        "_reduce_ground_in_drift(x::Number)",
                    "drift::Symbolics.Num" => "drift::Number",
                    "noise::Union{Nothing, Symbolics.Num}" => "noise::Union{Nothing, Number}",
                )
            end
            write(path, text)
        end
    end
    return src
end
