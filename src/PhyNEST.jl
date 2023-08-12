#Written by Sungsik Kong 2021-2022
#Last updated by Sungsik Kong 2023
module PhyNEST
    using DataFrames
    using DelimitedFiles
    using Distributed
    using CSV    
    #using ProgressMeter
    using Distributions
    using Optim, LineSearches
    using StatsBase
    using Dates
    using Logging
    using Suppressor
    import PhyloNetworks: 
        readTopology, 
        displayedTrees,
        Edge, 
        HybridNetwork,
        NNIRepeat!,
        PhyloNetworks.calculateNmov!,
        whichMove,
        readTopologyUpdate,
        writeTopologyLevel1,
        proposedTop!,
        printEverything,
        printEdges
    export 
        readTopology,
        displayedTrees,
        greet,
        checkDEBUG,
        Phylip,
        readPhylip, 
        show_sp,
        write_sp,
        storeCheckPoint,
        readCheckPoint,
        GetTrueProbsSymm,
        GetTrueProbsAsymm,
        GetTrueProbsAsymmTypes,
        sim_sp_counts,
        Network,
        quartets,
        get_quartets,
        printQuartets,
        get_start_theta,
        do_optimization,
        hill_climbing,
        simulated_annealing,
        phyne!,
        Dstat, showall

    include("check.jl")
    include("Objects.jl")
    include("readPhylip.jl")
    include("Probabilities.jl")
    include("Quartets.jl")
    include("MomentEstimat.jl")
    include("Theta.jl")
    include("Optimization.jl")
    include("Searches.jl")
    include("miscellaneous.jl")

end

#===
7.28 removed showwprogressmeter dependency; but looks a little less fancy
    modified some error messages

===#