# PhyNEST.jl

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sungsik-kong/PhyNEST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/stable)
[![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/dev)

## Overview

Julia package **[PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns** can:

- Read multi-locus sequence alignment (in a relaxed PHYLIP format) and parse site pattern frequencies for all permutations of four individuals (i.e., sequences).
- Estimate maximum composite likelihood species network using the parsed site pattern frequencies using hill climbing sor simulated annealing heuristic algorithms.

## Installation
To use PhyNEST, Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.

To install PhyNEST, run the following command in Julia prompt:
```@julia
using Pkg
Pkg.add("PhyNEST")
```
or
```@julia
] add PhyNEST
```
To check the installed version of PhyNEST, run the following command:
```@julia
] status PhyNEST
```

## A quick test
Once the installation is complete, load PhyNEST and execute function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@julia
using PhyNEST
greet()
```

## Documentation
Please see [here](https://sungsik-kong.github.io/PhyNEST.jl/) for detailed documentation.

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or make suggestions.


