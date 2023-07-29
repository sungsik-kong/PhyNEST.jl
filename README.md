# PhyNEST.jl

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/sungsik-kong/PhyNEST.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Stable-Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/stable)
[![Dev-Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://sungsik-kong.github.io/PhyNEST.jl/dev)

## Overview

Julia package **[PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns** can:

- Read multi-locus sequence alignment (in a relaxed PHYLIP format) and parse site pattern frequencies for all permutations of four individuals (i.e., sequences).
- Estimate maximum composite likelihood species network using the parsed site pattern frequencies using hill climbing or simulated annealing heuristic algorithms.

## Installation of the package PhyNEST
To use [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl), Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.

To install [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) run the following command in Julia prompt:
```@julia
julia> using Pkg
julia> Pkg.add("PhyNEST")
```
or
```@julia install
julia> ] add PhyNEST
```
Either command installs the package. To check the current version of installed [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl), run the following command:
```@julia
julia> ] status PhyNEST
```
## Test example
Once the installation is complete, load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) and run function `greet()` to see if the package has been added to the local machine. A greet message with current time should appear.
```@repl install
julia> using PhyNEST
julia> greet()
```

## Tutorials
A step-by-step [wiki tutorial](https://github.com/sungsik-kong/PhyNEST.jl/wiki) is available [here](https://github.com/sungsik-kong/PhyNEST.jl/wiki).

## Citation
Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on *BioRxiv* at [https://doi.org/10.1101/2022.11.14.516468](https://doi.org/10.1101/2022.11.14.516468).

## Documentation
Please see [here](https://sungsik-kong.github.io/PhyNEST.jl/) for detailed documentation.

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or make suggestions.