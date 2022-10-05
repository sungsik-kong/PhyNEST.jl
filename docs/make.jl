using Documenter
using PhyNE

makedocs(
    sitename = "PhyNE",
    format = Documenter.HTML(),
    modules = [PhyNE]
)

deploydocs(
    repo = "github.com/sungsik-kong/PhyNE.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
