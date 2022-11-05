# Installation

## Julia installation
Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single programming language. Current release of Julia is available [here](https://julialang.org/downloads/). Julia can be credited by citing the following paper:
>Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65–98. doi: [10.1137/141000671](10.1137/141000671).

## Installation of the package PhyNEST
To use [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl), Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.

To install [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) run the following command in Julia prompt:
```@julia
using Pkg
Pkg.add("PhyNEST")
```
or
```@julia install
] add PhyNEST
```
Either command installs the package. To check the current version of installed [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl), run the following command:
```@julia
] status PhyNEST
```
## Test example
Once the installation is complete, load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@repl install
using PhyNEST
greet()
```