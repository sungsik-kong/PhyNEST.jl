# PhyNE.jl

Documentation for PhyNE.jl: Estimating Maximum Composite Likelihood Phylogenetic Network

---

General info TBA


## References

TBA

# Manual
```@contents
Pages = [
    "manual/input.md",
    "manual/installation.md",
    "manual/quartet.md",
    "manual/optimization.md",
    "manual/networkest.md",
    "manual/postanalysis.md",
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