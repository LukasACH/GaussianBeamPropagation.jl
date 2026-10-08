using Documenter, GaussianBeamPropagation

makedocs(;
    sitename="GaussianBeamPropagation",
    pages=Any[
        "Introduction"=>"index.md",
        "How to get started"=>"getting_started.md",
        "API"=>Any[
            "Helpers"=>"api/helpers.md",
            "Elements"=>"api/elements.md",
        ],
        "Integration with other packages"=>"integration.md",
        "API Reference"=>"reference.md",
    ],
)

deploydocs(; repo="github.com/LukasACH/GaussianBeamPropagation.jl.git")
