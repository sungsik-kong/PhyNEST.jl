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
    get_quartets
    get_leaf_name
    get_leaf_number
    list_all_quartets
    GetChild
    GetParent
    childnode
    parentnode
    mrca
    get_most_recent_common_ancestors
    identify_sym_type
    tauNum
    backtauNum
    get_unique_tau_labels
    dictionary_phylip
    dictionary_topology
    dictionary_combined
    get_matching_quartet
    transfer_site_pattern_frequencies
    get_gamma
    printQuartets
    GetTrueProbsSymm
    GetTrueProbsAsymm
    GetTrueProbsAsymmTypes
    GetTrueProbsNetTypes
    sim_sp_counts
    do_optimization
    get_start_gamma_info
    get_parent_child
    get_initial_taus
    backTransform
    backtransform_parameters
    update_topology
    hill_climbing
    burn_in
    simulated_annealing
    initiate_search
    phyne!
    get_average
    greet
    checkDEBUG
    Dstat
    HyDe
```