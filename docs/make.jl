using Documenter
#using PhyNEST

makedocs(
	sitename = "PhyNEST",
	pages = [
		"index.md",
		"Manual" => [
			#"manual/installation.md",
			#"manual/input.md",
			#"manual/networkest.md",
#			"manual/quartet.md"
#			"manual/others.md"
		]
	]
)

deploydocs(
	repo = "github.com/sungsik-kong/PhyNEST.jl.git",
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
