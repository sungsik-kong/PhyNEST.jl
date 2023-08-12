# PhyNEST

Documentation for a Julia package **[PhyNEST](https://github.com/sungsik-kong/PhyNEST.jl)**: Phylogenetic Network Estimation using SiTe patterns.

---
## Tutorials
A step-by-step [wiki tutorial](https://github.com/sungsik-kong/PhyNEST.jl/wiki) is available, that has been done for the 2023 Botany workshop.

## Reference
Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on *BioRxiv* at [https://doi.org/10.1101/2022.11.14.516468](https://doi.org/10.1101/2022.11.14.516468).

## Getting help
Please use [google group](https://groups.google.com/g/phynest-users) to report bugs or post questions and/or suggestions.


## Functions
```@docs
    PhyNEST.greet
    PhyNEST.readPhylip
    PhyNEST.readPhylipFile
    PhyNEST.PhylipFileInfo
    PhyNEST.getUniqueQuartets
    PhyNEST.sitePatternCounts
    PhyNEST.spRearrange
    PhyNEST.show_sp
    PhyNEST.sitePatternsToDF
    PhyNEST.write_sp
    PhyNEST.storeCheckPoint
    PhyNEST.readCheckPoint

    PhyNEST.Dstat
    PhyNEST.showall
```
## Types
```@docs
    Phylip
```

## Index
```@index
    
```