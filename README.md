# PhyNEST: Phylogenetic Network Estimation using SiTe patterns

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sungsik-kong/PhyNEST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/dev)

## Overview

Julia package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) can:

- Read multi-locus sequence alignment (in a relaxed Phylip format) and parse site pattern frequencies for all possible permutations of four individuals (i.e., quartet).
- Estimate maximum composite likelihood species networks from the parsed site patterns extracted using hill climbing or simulated annealing heuristic algorithms.

## Installation
To install PhyNEST, run the following command in Julia prompt:
```@julia install
using Pkg
Pkg.add("PhyNEST")
```
or
```@julia install
] add PhyNEST
```
Either command installs the package.

To check the installed version of PhyNEST, run the following command:
```@julia install
] status PhyNEST
```
## A quick test
Once the installation is complete, execute PhyNEST and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@repl install
using PhyNEST
greet()
```

## Documentation
Please see [here](https://sungsik-kong.github.io/PhyNEST.jl/) for detailed documentation.

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or make suggestions.


