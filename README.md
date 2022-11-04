# PhyNEST.jl

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sungsik-kong/PhyNEST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/stable)
[![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/dev)

## Overview

Julia package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns can:

- Read multi-locus sequence alignment (in a relaxed Phylip format) and parse site pattern frequencies for all possible permutations of four individuals (i.e., quartet).
- Estimate maximum composite likelihood species networks from the parsed site patterns extracted using hill climbing or simulated annealing heuristic algorithms.

## Installation
To install PhyNEST, run the following command in Julia prompt:
```@julia
using Pkg
Pkg.add("PhyNEST")
```
or
```@julia
] add PhyNEST
```
Either command installs the package.

To check the currently installed version of PhyNEST, run the following command:
```@julia
] status PhyNEST
```
## A quick test
Once the installation is complete, execute PhyNEST and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@julia
using PhyNEST
greet()
```

## Documentation
Please see [here](https://sungsik-kong.github.io/PhyNEST.jl/) for detailed documentation.

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or make suggestions.


