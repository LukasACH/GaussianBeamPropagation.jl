using Documenter, GaussianBeamPropagation

makedocs(;
    sitename="GaussianBeamPropagation",
    pages=Any[
        "Home"=>"index.md",
    ],
)
deploydocs(
    ; repo="github.com/LukasACH/GaussianBeamPropagation.jl.git",
)
