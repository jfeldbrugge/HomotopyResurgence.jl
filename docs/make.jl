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
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/HomotopyResurgence.jl",
    devbranch="main",
)
