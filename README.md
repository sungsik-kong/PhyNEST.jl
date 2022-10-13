# PhyNE: Estimating Maximum Composite Likelihod Phylogenetic Network

[![Build Status](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/sungsik-kong/PhyNE.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Overview

Documentation for a julia package [PhyNE.jl](https://github.com/sungsik-kong/PhyNE.jl): **Phy**logenetic **N**etwork **E**stimation using composite likelihood directly from the sequence data using site pattern frequency. 

PhyNe (**Phy**logenetic **N**etwork **e**stimation) is a Julia package that can:

- Read multilocus sequence alignment (in relaxed phylip format) and extract site pattern frequencies for all possible permutations of quartet (i.e., four individuals).
- Extract rooted quartet trees from phylogenetic trees and networks in newick (or extended newick for networks) format.
- Compute expected site pattern frequencies for a quartet tree given the speciation times, effective population size parameter, and mutation rate.
- Identify three speciation times for a quartet given the site pattern frequencies and effective population size parameter.
- Estimate species networks directly from multilocus sequence alignment using hill-climbing or simulated annealing algorithms.
- Detect hybridization using HyDe and the D-statistic.

## Quick Tutorial

## Getting help
Please use [google group](https://groups.google.com/g/phyne-users) to report bugs or make suggestions.
