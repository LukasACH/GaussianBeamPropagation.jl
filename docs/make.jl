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
        "Fitting"=>"fitting.md",
        "Integration with other packages"=>"integration.md",
        "API Reference"=>"reference.md",
        "Future work"=>Any[
            "Optical assemblies"=>"future_work/optical_assembly.md",
            "Propagation along a system"=>"future_work/propagation.md",
            "Moment types"=>"future_work/moment_types.md",
        ],
    ],
)

deploydocs(; repo="github.com/LukasACH/GaussianBeamPropagation.jl.git", devbranch="main")
