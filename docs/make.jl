using Homotopy
using Documenter

DocMeta.setdocmeta!(Homotopy, :DocTestSetup, :(using Homotopy); recursive=true)

makedocs(;
    modules=[Homotopy],
    authors="Job Feldbrugge <14946916+jfeldbrugge@users.noreply.github.com> and contributors",
    sitename="Homotopy.jl",
    format=Documenter.HTML(;
        canonical="https://jfeldbrugge.github.io/Homotopy.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jfeldbrugge/Homotopy.jl",
    devbranch="main",
)
