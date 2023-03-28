# PhyNEST.jl

Documentation for a Julia package **[PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns**.

---

## References
Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on *BioRxiv* at [https://doi.org/10.1101/2022.11.14.516468](https://doi.org/10.1101/2022.11.14.516468).

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or post questions and/or suggestions.

## Manual
```@contents
Pages = [
    "manual/installation.md",
    "manual/input.md",
    "manual/networkest.md",
    ]
```

## Functions
```@docs
    greet
    readPhylip
    show_sp
    write_sp
    storeCheckPoint
    readCheckPoint
    GetTrueProbsSymm
    GetTrueProbsAsymm
    GetTrueProbsAsymmTypes
    sim_sp_counts
    get_quartets
    printQuartets
    get_start_theta
    do_optimization
    hill_climbing
    simulated_annealing
    phyne!
```
## Types
```@docs
    Phylip
    Network
    quartets
```

## Index
```@index
    readPhylip
    readPhylipFile
    PhylipFileInfo
    getUniqueQuartets
    sitePatternCounts
    spRearrange
    show_sp
    sitePatternsToDF
    write_sp
    writeSitePatternCounts
    storeCheckPoint
    readCheckPoint
```