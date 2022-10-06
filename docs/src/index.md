# PhyNE.jl

Documentation for a julia package [PhyNE.jl](https://github.com/sungsik-kong/PhyNE.jl): **Phy**logenetic **N**etwork **E**stimation using composite likelihood directly from the sequence data using site pattern frequency. 

---


## References

TBA

## Getting help
Please use [google group](https://groups.google.com/g/phyne-users) to report bugs or make suggestions.

## Manual
```@contents
Pages = [
    "manual/installation.md",
    "manual/input.md",
    "manual/quartet.md",
    "manual/networkest.md",
    "manual/others.md",
    ]
```

## Functions

```@docs
    readPhylipFile!
    readCheckPoint
    writeSitePatternCounts
    mrca
    nature
    extractNQuartets
    printQuarts
    GetTrueProbsSymm
    GetTrueProbsAsymm
    GetTrueProbsNetTypes
    simspcounts
    momentEstimat
    startTheta
    Optimization
    hillClimb
    maxburnin
    simAnneal

    PhyNe!

    Dstat
    HyDe
    LRT
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
    printQuarts
    Phylip
    mrca
    nature
    nquartets
    Nquartets
    PhyNe!
```