var documenterSearchIndex = {"docs":
[{"location":"#PhyNEST","page":"PhyNEST","title":"PhyNEST","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Documentation for a Julia package PhyNEST: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"","category":"page"},{"location":"#Tutorials","page":"PhyNEST","title":"Tutorials","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"A step-by-step wiki tutorial is available, that has been done for the 2023 Botany workshop.","category":"page"},{"location":"#Reference","page":"PhyNEST","title":"Reference","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.","category":"page"},{"location":"#Getting-help","page":"PhyNEST","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Functions","page":"PhyNEST","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    PhyNEST.greet\n    PhyNEST.readPhylip\n    PhyNEST.readPhylipFile\n    PhyNEST.PhylipFileInfo\n    PhyNEST.getUniqueQuartets\n    PhyNEST.sitePatternCounts\n    PhyNEST.spRearrange\n    PhyNEST.show_sp\n    PhyNEST.sitePatternsToDF\n    PhyNEST.write_sp\n    PhyNEST.storeCheckPoint\n    PhyNEST.readCheckPoint\n\n    PhyNEST.GetTrueProbsSymm\n    PhyNEST.GetTrueProbsAsymm\n    PhyNEST.GetTrueProbsAsymmTypes\n    PhyNEST.sim_sp_freq\n\n    PhyNEST.phyne!\n\n    PhyNEST.Dstat\n    PhyNEST.showallDF","category":"page"},{"location":"#PhyNEST.greet","page":"PhyNEST","title":"PhyNEST.greet","text":"greet()\n\nDisplays a greeting with citation information. No argument needed. \n\nExample\n\njulia> greet()\nThank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.\nPlease report bugs or make suggestions to https://groups.google.com/g/phynest-users.\nIf you conduct an analysis using PhyNEST, please cite:\nSungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood.\nPreprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readPhylip","page":"PhyNEST","title":"PhyNEST.readPhylip","text":"readPhylip(inputfile::AbstractString; \n            writecsv=false::Bool,\n            csvname=\"\"::AbstractString,\n            showProgress=true::Bool,\n            tdigits=3::Integer,\n            checkpoint=false::Bool)\n\nImport, read, and parse the input phylip file. File name must be specified as string. By default, progress bar is displayed, and .csv and .ckp files are NOT produced.  See optional arguments below and modify if needed.\n\nMandatory argument\n\ninputfile     Name of the phylip file as a string\n\nOptional arguments\n\nwritecsv       (default=false) A boolean arguement to allow writing site pattern frequencies in .csv file\ncsvname        (default=sitePatternCounts_inputfile.csv) A string that will be name of the .csv file \nshowProgress   (default=true) A boolean argument for visualizing the progress of alignment parsing\ntdigits        (default=3) Number of decimal points to record timme taken to parse the file in seconds\ncheckpoint     (default=false) A boolean argument to store the Phylip object as a .ckp file. (Warning: .ckp file can be large)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readPhylipFile","page":"PhyNEST","title":"PhyNEST.readPhylipFile","text":"readPhylipFile(inputfile::AbstractString,\n                writecsv::Bool,\n                csvname::AbstractString,\n                showProgress::Bool)\n\nExceuted while running the function readPhylip.  Fills in the attributes in the Phylip object, except for the Phylip.time. \n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.PhylipFileInfo","page":"PhyNEST","title":"PhyNEST.PhylipFileInfo","text":"PhylipFileInfo(inputfile::AbstractString, \n                p::Phylip, \n                showProgress::Bool)\n\nFunction that shortens?summarizes? the sequence alignment into two matrices: \n\nThe one with unique sites and \nAnother with how many times each column in the previous matrix occurs throughout the alignment.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.getUniqueQuartets","page":"PhyNEST","title":"PhyNEST.getUniqueQuartets","text":"getUniqueQuartets(p::Phylip)\n\nGetting some combinations for four sequences (quartets). This will result in n! quartets where n is the number of sequences. I believe a better way to do this exists.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sitePatternCounts","page":"PhyNEST","title":"PhyNEST.sitePatternCounts","text":"sitePatternCounts(p::Phylip,\n                ppbase::Array,\n                counts::Array)\n\nComputes observed site pattern frequencies from UniqueBase and BaseCounts obtained from the function PhylipFileInfo.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.spRearrange","page":"PhyNEST","title":"PhyNEST.spRearrange","text":"spRearrange(p::Phylip)\n\nRearranges the elemetns of each unique quartet obtained using getUniqueQuartets, and also suffles the site pattern frequencies accordingly. This prevents 24 redundant computations, saving computation time.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.show_sp","page":"PhyNEST","title":"PhyNEST.show_sp","text":"show_sp(p::Phylip)\n\nPretty name for displaying all quartet and the corresponding site pattern frequencies on screen in the DataFrame format. It may result in excessively long table when there are many sequences in the input file.\n\nMandatory argument\n\np     Phylip object\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sitePatternsToDF","page":"PhyNEST","title":"PhyNEST.sitePatternsToDF","text":"sitePatternsToDF(p::Phylip)\n\nExtracts the quartet and site pattern information from the Phyliip object, and reorganizes them in the DataFrame format.  This function does all hard work work for show_sp.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.write_sp","page":"PhyNEST","title":"PhyNEST.write_sp","text":"write_sp(p::Phylip; \n        csvname=\"PhyNEST_sp\"::AbstractString)\n\nPretty name for exporting all quartet and the corresponding site pattern frequencies in a .csv format.\n\nMandatory argument\n\np         Phylip object \n\nOptional argument\n\ncsvname   Filename for the .csv output can be given, otherwise we use PhyNEST_sp.csv by default.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.storeCheckPoint","page":"PhyNEST","title":"PhyNEST.storeCheckPoint","text":"storeCheckPoint(p::Phylip)\n\nStores the phylip object into a .ckp file so a user does not have to repeat the phylip parsing again if working with the  same dataset next time. By default it creates a .ckp file that has the same filename as the input phylip file. File size can get pretty large... \n\nMandatory argument\n\np         Phylip object\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readCheckPoint","page":"PhyNEST","title":"PhyNEST.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file and creates a filled in phylip object.\n\nMandatory argument\n\nckpfile   Name of the checkpoint file\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsSymm","page":"PhyNEST","title":"PhyNEST.GetTrueProbsSymm","text":"GetTrueProbsSymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64)\nGetTrueProbsSymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64,\n                alpha::Float64)\n\nComputes true site pattern probabilities for the symmetric quartet tree: ((1,2),(3,4));. The fifteen quartet site pattern probabilities are returned  in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nThree speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified.  See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nMandatory arguments\n\nmyt1      Speciation time for the common ancestor of species 1 and 2 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\nOptional argument\n\nalpha (dafault=4/3)\n\nExample\n\njulia> symprob=GetTrueProbsSymm(1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.2404720290824025\n 0.00047873583861244776\n 0.00047873583861244776\n 0.0007310388098746822\n 3.639464668511165e-6\n 0.0007211284736310421\n 6.790832193120881e-6\n 1.461658030751117e-6\n 6.790832193120881e-6\n 0.0007211284736310421\n 1.461658030751117e-6\n 1.461658030751117e-6\n 1.461658030751117e-6\n 6.304951604760085e-6\n 2.955516269578535e-8\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsAsymm","page":"PhyNEST","title":"PhyNEST.GetTrueProbsAsymm","text":"GetTrueProbsAsymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64)\nGetTrueProbsAsymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64,\n                alpha::Float64)\n\nComputes true site pattern probabilities for the asymmetric quartet tree: (1,(2,(3,4)));. The fifteen quartet site pattern probabilities are returned  in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nThree speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified.  See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nMandatory arguments\n\nmyt1      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 2 and (3,4) in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\n##Optional argument\n\nalpha (dafault=4/3)\n\nExample\n\njulia> asymprob=GetTrueProbsAsymm(1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.24055000044773364\n 0.00045274538350207147\n 0.00045274538350207147\n 0.00027457470635933667\n 1.7612660121311035e-6\n 0.0006951380185206657\n 3.283594866865977e-5\n 1.4343273481706927e-6\n 3.283594866865977e-5\n 0.0011794138135323934\n 2.401907811548733e-6\n 1.4343273481706927e-6\n 2.401907811548733e-6\n 5.394333411761533e-6\n 2.7254257479319977e-8\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsAsymmTypes","page":"PhyNEST","title":"PhyNEST.GetTrueProbsAsymmTypes","text":"GetTrueProbsAsymmTypes(type::Integer,\n                        myt1::Float64,\n                        myt2::Float64,\n                        myt3::Float64,\n                        theta::Float64)\nGetTrueProbsAsymmTypes(type::Integer,\n                        myt1::Float64,\n                        myt2::Float64,\n                        myt3::Float64,\n                        theta::Float64,\n                        alpha::Float64)\n\nComputes true site pattern probabilities for any of the four the asymmetric quartet trees: \n\nType 1: (1,((2,3),4));\nType 2: ((1,(2,3)),4); \nType 3: (1,(2,(3,4))); \nType 4: (((1,2),3),4);\n\nThe fifteen quartet site pattern probabilities are returned in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nThree speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified.  See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nMandatory arguments\n\ntype      Type of the asymmetric quartet in interest. See above.\nmyt1      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 2 and (3,4) in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\nOptional argument\n\nalpha (dafault=4/3)\n\nExample\n\njulia> GetTrueProbsAsymmTypes(1,1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.24055000044773364\n 0.0006951380185206657\n 0.00045274538350207147\n 3.283594866865977e-5\n 1.4343273481706927e-6\n 0.00045274538350207147\n 3.283594866865977e-5\n 1.4343273481706927e-6\n 0.00027457470635933667\n 0.0011794138135323934\n 5.394333411761533e-6\n 1.7612660121311035e-6\n 2.401907811548733e-6\n 2.401907811548733e-6\n 2.7254257479319977e-8\n\njulia> GetTrueProbsAsymmTypes(2,1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.24055000044773364\n 0.0011794138135323934\n 0.00045274538350207147\n 3.283594866865977e-5\n 2.401907811548733e-6\n 0.00045274538350207147\n 3.283594866865977e-5\n 2.401907811548733e-6\n 0.00027457470635933667\n 0.0006951380185206657\n 5.394333411761533e-6\n 1.7612660121311035e-6\n 1.4343273481706927e-6\n 1.4343273481706927e-6\n 2.7254257479319977e-8\n\njulia> GetTrueProbsAsymmTypes(3,1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.24055000044773364\n 0.00045274538350207147\n 0.00045274538350207147\n 0.00027457470635933667\n 1.7612660121311035e-6\n 0.0006951380185206657\n 3.283594866865977e-5\n 1.4343273481706927e-6\n 3.283594866865977e-5\n 0.0011794138135323934\n 2.401907811548733e-6\n 1.4343273481706927e-6\n 2.401907811548733e-6\n 5.394333411761533e-6\n 2.7254257479319977e-8\n\njulia> GetTrueProbsAsymmTypes(4,1.0,2.0,3.0,0.003,4/3)\n15-element Vector{Float64}:\n 0.24055000044773364\n 0.0011794138135323934\n 0.0006951380185206657\n 0.00027457470635933667\n 5.394333411761533e-6\n 0.00045274538350207147\n 3.283594866865977e-5\n 2.401907811548733e-6\n 3.283594866865977e-5\n 0.00045274538350207147\n 2.401907811548733e-6\n 1.4343273481706927e-6\n 1.4343273481706927e-6\n 1.7612660121311035e-6\n 2.7254257479319977e-8\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sim_sp_freq","page":"PhyNEST","title":"PhyNEST.sim_sp_freq","text":"sim_sp_freq(type::Integer,\n            myt1::Float64,\n            myt2::Float64,\n            myt3::Float64,\n            theta::Float64)\nsim_sp_freq(type::Integer,\n            myt1::Float64,\n            myt2::Float64,\n            myt3::Float64,\n            theta::Float64,\n            alpha::Float64,\n            length::Integer)\n\nGenerates site pattern frequencies for the five possible quartet topologies modeled as  a multinomial random variable under the assumption that the observed sites are independent,  conditional on the species tree. The five topologies are:\n\nType 0: ((1,2),(3,4));\nType 1: (1,((2,3),4));\nType 2: ((1,(2,3)),4); \nType 3: (1,(2,(3,4))); \nType 4: (((1,2),3),4).\n\nThe fifteen quartet site pattern frequencies are returned in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nMandatory arguments\n\ntype      Type of the asymmetric quartet in interest. See above.\nmyt1      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 2 and (3,4) in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\n##Optional argument\n\nalpha (dafault=4/3)\nlength (default=1000000) Alignment length\n\n##Example\n\njulia> sim_sp_freq(1,1.0,2.0,3.0,0.003,4/3,50000)\n15-element Vector{Int64}:\n 48073\n   411\n   261\n    14\n     6\n   285\n    33\n     1\n   171\n   729\n     9\n     1\n     2\n     4\n     0\njulia> sim_sp_freq(2,1.0,2.0,3.0,0.003)\n15-element Vector{Int64}:\n 962053\n  14160\n   5460\n    410\n     58\n   5405\n    380\n     57\n   3372\n   8434\n    122\n     32\n     32\n     24\n      1     \n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.phyne!","page":"PhyNEST","title":"PhyNEST.phyne!","text":"phyne!(starting_topology::HybridNetwork,p::Phylip,outgroup::String;\n        hmax=1::Integer,\n        maximum_number_of_steps=250000::Integer,\n        maximum_number_of_failures=100::Integer,\n        number_of_itera=1000::Integer,\n        number_of_runs=10::Integer,\n        do_hill_climbing=true::Bool,\n        number_of_burn_in=25::Integer,\n        k=10::Integer,\n        cons=0.5::Float64,\n        alph=0.5::Float64,\n        filename=\"\"::AbstractString)\n\nphyne! is fine. phyne! executes function initiate_search(args).\n\nEstimate the species network (or tree if hmax=0) using maximum composite likelihood. The search begins from the starting_topology which can be either estimated or randomly generated. Starting topology can be either tree or a network with <=hamx. Outgroup taxon  must be specified to root the network. \n\nMandatory arguments\n\nstarting_topology       Starting topology in HybridNetwork object created using the function readTopology \np                       Sequence alignment parsed in Phylip object using the function readPhylip\noutgroup                Name of the outgroup taxa\n\nOptional arguments\n\nGeneric\n\nhmax (dafault=1)                            Maximum number of hybridizations to be included in the final network\ndo_hill_climbing (default=true)             When =true, network is searched using hill climbing and when =false, it searches using simulated annealing\nnumber_of_runs (default=10)                 Number of independent runs\n\nNetwork space search\n\nmaximum_number_of_steps (default=250000)    Maximum number of steps before the search terminates\nmaximum_number_of_failures (default=100)    Maximum number of consecutive failures (i.e., rejecting the proposed topology) before the search teminates\n\nOptimization\n\nnumber_of_itera (default=1000)\n\nFor simulated annealing\n\nnumber_of_burn_in (default=25)              \nk (default=10)                              Specifies the number of best networks to be stored at each run\ncons (default=0.5)                          \nalph (default=0.5)\n\nOutput\n\nfilename (default=\"\") Specifies the name of the output file. If unspecified, it will use PhyNEST_hc or PhyNEST_sa depending on the heuristic method applied.\nThe best network estimated throughout the entire runs is written in the extended Newick format and stored in file with an extension .out. \n\nFirst line in the file is readable by PhyloNetworks that can be visualized using PhyloPlots. Third (last) line is the identical network but readable by DendroScope. \n\nSpecifics of the network search from each independent run is stored in the file with an extension .log. \n\nIt summarizes how many of each network 'moves' were made during the search, the reason for terminating the search  (e.g., researched maximum number of steps, reached maximum number of consecutive failures, or else), and  the final network estimated and its composite likelihood. When do_hill_climbing=true, a single best network is selected for each run, but when do_hill_climbing=false it prints k best networks visited.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.Dstat","page":"PhyNEST","title":"PhyNEST.Dstat","text":"Dstat(outgroup::String, p::Phylip;\n    pval=0.05::Float64, \n    displayall=false::Bool,\n    writecsv=false::Bool, \n    filename=\"\"::AbstractString)\n\nConducts Patterson's D-statistic test. The result prints site pattern frequencies ABAB and ABBA used to compute  the D-statistic, Z-score, and the p-value for each quartet tested. Significance is marked with an asterisk. Function showall(df) can be subsequently used to show all rows.\n\nMandatory arguments\n\noutgroup     Name of the outgroup taxa\np   The Phylip object\n\nOptional arguments\n\npval       (default=0.05) Alpha level for significance\ndisplay_all (default=false) If set as true, the function print test results for every quartet. By default, it only prints those quartets where signficance was found.\nwritecsv (default=false) If true, the result is stored in .csv file in the working directory\nfilename Specifies .csv file name if writecsv=true. If unspecified, the result is stored as Dstat-out.csv\n\nExample\n\njulia> p=readPhylip(\"sample_n4h1.txt\")\njulia> df=Dstat(\"4\",p)\nTip: if neccessary, use showall(df) function to see all the rows.\n2×10 DataFrame\n Row │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat     Zscore   pvalue   significance\n     │ String    String  String  String  Int64  Int64  Float64   Float64  Float64  String\n─────┼──────────────────────────────────────────────────────────────────────────────────────────\n   1 │ 4         3       2       1        1427   7852  0.692424  66.6995      0.0  *\n   2 │ 4         1       2       3        1427   7836  0.691892  66.5908      0.0  *\n\njulia> df=Dstat(\"4\",p,display_all=true)\nTip: if neccessary, use showall(df) function to see all the rows.\n6×10 DataFrame\nRow │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat        Zscore      pvalue    significance\n    │ String    String  String  String  Int64  Int64  Float64      Float64     Float64   String\n────┼─────────────────────────────────────────────────────────────────────────────────────────────────\n  1 │ 4         3       1       2        7852   1427  -0.692424    -66.6995    1.0\n  2 │ 4         3       2       1        1427   7852   0.692424     66.6995    0.0       *\n  3 │ 4         1       3       2        7836   1427  -0.691892    -66.5908    1.0\n  4 │ 4         1       2       3        1427   7836   0.691892     66.5908    0.0       *\n  5 │ 4         2       3       1        7836   7852   0.00101989    0.127743  0.449176\n  6 │ 4         2       1       3        7852   7836  -0.00101989   -0.127743  0.550824\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.showallDF","page":"PhyNEST","title":"PhyNEST.showallDF","text":"showallDF(df::DataFrame)\n\nPrint all rows of the DataFrame object using the package CSV.    \n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNEST","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    Phylip","category":"page"},{"location":"#PhyNEST.Phylip","page":"PhyNEST","title":"PhyNEST.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type with the following attributes:\n\nfilename      Name of the input phylip alignment file\n\ntime          Time taken for parsing the input in seconds\n\nnumtaxa       Number of taxa given at the first line of the input file\n\nseqleng       Sequence length given at the first line of the input file\n\nnametaxa      Sequence names given in the input file\n\ncounttaxa     Unique integer identifier given to each individual in the order of appearance in the input file\n\nallquartet    All combinations of quartets for counttaxa\n\nindex         Just convert allquartet elements into arbitrary index numbers\n\nspcounts      Arrays of 15 site pattern frequencies for each quartet in allquartet\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNEST","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    ","category":"page"}]
}
