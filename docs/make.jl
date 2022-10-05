using Documenter
using PhyNE

makedocs(
    sitename = "PhyNE",
    format = Documenter.HTML(),
    modules = [PhyNE]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
