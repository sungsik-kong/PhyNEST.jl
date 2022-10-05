# PhyNE.jl

Documentation for a julia package [PhyNE.jl](https://github.com/sungsik-kong/PhyNE.jl): Estimating Maximum Composite Likelihood Phylogenetic Network directly from the sequence data using site pattern frequency. 

---

Please report bugs or suggestions to the [google group](https://groups.google.com/g/phyne-users). Some frequently asked questions may are already answered in the group.


## References

TBA

# Manual
```@contents
Pages = [
    "manual/input.md",
    "manual/installation.md",
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