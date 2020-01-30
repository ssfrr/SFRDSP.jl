using Documenter, SFRDSP

makedocs(;
    modules=[SFRDSP],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/ssfrr/SFRDSP.jl/blob/{commit}{path}#L{line}",
    sitename="SFRDSP.jl",
    authors="Spencer Russell",
    assets=String[],
)

deploydocs(;
    repo="github.com/ssfrr/SFRDSP.jl",
)
