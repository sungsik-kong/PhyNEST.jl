# PhyNEST.jl

Documentation for a julia package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl): **Phy**logenetic **N**etwork **E**stimation using *S*i*T*e patterns.

---

## References
TBA

## Getting help
Please use [google group](https://groups.google.com/g/phyne-users) to report bugs or post questions and/or suggestions.

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
    readPhylipFile!
    readCheckPoint
    writeSitePatternCounts
    extractNQuartets
    extractQuartets
    printQuartets
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
    Phylip
    nquartets
    Nquartets
    PhyNE!
```
