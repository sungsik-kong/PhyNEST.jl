# PhyNEST.jl

Documentation for a Julia package **[PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): Phylogenetic Network Estimation using SiTe patterns**.

---

## References
TBA

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or post questions and/or suggestions.

## Manual
```@contents
Pages = [
    "manual/installation.md",
    "manual/input.md",
    "manual/networkest.md",
    "manual/quartet.md"
    ]
```
## Functions
```@docs
    readPhylipFile!
    readCheckPoint
    writeSitePatternCounts
    extractNQuartets
    extractQuartets
    printQuartets
    TrueSitePatternSymm
    TrueSitePatternAsymm
    simspcounts
    PhyNE!
```
## Types
```@docs
    Phylip
    nquartets
    Nquartets
```

## Index
```@index
    readPhylipFile!
    readCheckPoint
    writeSitePatternCounts
    extractNQuartets
    printQuartets
    TrueSitePatternSymm
    TrueSitePatternAsymm
    PhyNE!
    Phylip
    nquartets
    Nquartets
```
