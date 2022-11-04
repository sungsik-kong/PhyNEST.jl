# Installation

## Julia installation

Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single language. Current release of Julia is available [here](https://julialang.org/downloads/). Julia can be creditted by citing the following paper:

Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: 10.1137/141000671.

## Installation of the package PhyNE
To install the package, type in Julia prompt:
```@julia install
using Pkg
Pkg.add("PhyNEST")
``

## Test example
Once the installation is complete, execute PhyNEST and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@repl install
using PhyNEST
greet()
``