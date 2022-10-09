using Documenter
using PhyNE

makedocs(
	sitename = "PhyNE.jl",
	pages = [
		"index.md",
		"Manual" => [
			"manual/installation.md",
			"manual/input.md",
			"manual/quartet.md",
			"manual/networkest.md",
			"manual/others.md"
		]
	]
)

deploydocs(
	repo = "github.com/sungsik-kong/PhyNE.jl.git",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
