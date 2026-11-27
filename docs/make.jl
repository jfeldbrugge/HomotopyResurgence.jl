ENV["JULIA_PROGRESS"] = "false"

using HomotopyResurgence
using Documenter

DocMeta.setdocmeta!(HomotopyResurgence, :DocTestSetup, :(using HomotopyResurgence); recursive=true)

makedocs(;
    modules=[HomotopyResurgence],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="HomotopyResurgence.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/HomotopyResurgence.jl",
        edit_link="main",
        assets=String[],
        size_threshold=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Theory" => "theory.md",
        "Tutorial" => "tutorial.md",
        "Citations" => "citations.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/HomotopyResurgence.jl",
    push_preview = true,
    devbranch="main",
)
