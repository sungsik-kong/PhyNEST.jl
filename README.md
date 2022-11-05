# PhyNEST.jl

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sungsik-kong/PhyNEST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/stable)
[![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/dev)

## Overview

Julia package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns can:

- Read multi-locus sequence alignment (in a relaxed `PHYLIP` format) and parse site pattern frequencies for all permutations of four individuals (i.e., sequences).
- Estimate maximum composite likelihood species networks using the parsed site pattern frequencies using hill climbing sor simulated annealing heuristic algorithms.

## Installation of the package PhyNEST
To install [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) run the following command in Julia prompt:
```@julia install
using Pkg
Pkg.add("PhyNEST")
```
or
```@julia install
] add PhyNEST
```
Either command installs the package.

To check the current version of installed [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl), run the following command:
```@julia install
] status PhyNEST
```
## Test example
Once the installation is complete, load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@repl install
using PhyNEST
greet()
```

## Documentation
Please see [here](https://sungsik-kong.github.io/PhyNEST.jl/) for detailed documentation.

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or make suggestions.