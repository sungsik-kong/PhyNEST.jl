# PhyNEST

Documentation for a Julia package **[PhyNEST](https://github.com/sungsik-kong/PhyNEST.jl)**: Phylogenetic Network Estimation using SiTe patterns.

---
## Tutorials
A step-by-step [wiki tutorial](https://github.com/sungsik-kong/PhyNEST.jl/wiki) for the major functions in `PhyNEST` is available.

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

    PhyNEST.GetTrueProbsSymm
    PhyNEST.GetTrueProbsAsymm
    PhyNEST.GetTrueProbsAsymmTypes
    PhyNEST.sim_sp_freq

    PhyNEST.phyne!

    PhyNEST.Dstat
    PhyNEST.showallDF
    PhyNEST.HyDe
```
## Types
```@docs
    Phylip
```

## Index
```@index
    
```