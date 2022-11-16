var documenterSearchIndex = {"docs":
[{"location":"manual/input/#Input","page":"Input","title":"Input","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Current version of PhyNEST.jl takes a relaxed, sequential, PHYLIP formatted DNA alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"A sample DNA alignment file n5h1_5k.txt is provided in the /example folder of the package PhyNEST.jl (or can be downloaded here). The alignment contains five sequences 1, 2, 3, 4, and 5 and is 2,500,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): (5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/input/#Parsing-DNA-alignment-data","page":"Input","title":"Parsing DNA alignment data","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Load PhyNEST.jl in Julia and move to the working directory that contains the alignment. The alignment is parsed with the function readPhylipFile!(). The progress can be visualized with showProgress option that takes a boolean argument (set as false in this example for brevity). ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Once the parsing the alignment is complete, a checkpoint file n5h1_5k.txt.ckp should appear in the directory. Note the .ckp file will have the same name as the alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"using PhyNEST\ndatapath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"n5h1_5k.txt\");\nphydata = readPhylipFile!(datapath, showProgress=false)","category":"page"},{"location":"manual/input/#Reading-.ckp-file","page":"Input","title":"Reading .ckp file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"In case of multiple network analyses using the same alignment, a user can bypass the data parsing process by importing the .ckp file that is created after parsing the alignment for the first time. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Reading a checkpoint file using the function readCheckPoint() as shown below.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"ckppath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"n5h1_5k.txt.ckp\");\nckpdata = readCheckPoint(ckppath)","category":"page"},{"location":"manual/input/#Exporting-as-.csv-file","page":"Input","title":"Exporting as .csv file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"The observed quartet site pattern frequencies can be exported in .csv file. When parsing the alignment using the function readPhylipFile!() set the option writecsv=true and the .csv file is created upon the completion of the parsing. The name of the .csv file can be set using the argument csvname. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Here, we set the name as tutorial_n5h1. The produced .csv file can be visualized using function readCSVFile(csvfile). ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"phydata = readPhylipFile!(datapath, showProgress=false,\n                        writecsv=true, csvname=\"tutorial_n5h1\")\ncsvpath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"tutorial_n5h1.csv\");\ncsvdf = readCSVFile(csvpath)","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"If a user forgot to set writecsv=true, the function writeSitePatternCounts can be used to export. See here for more information.","category":"page"},{"location":"manual/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual/installation/#Julia-installation","page":"Installation","title":"Julia installation","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single programming language. Current release of Julia is available here. Julia can be credited by citing the following paper:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65–98. doi: 10.1137/141000671.","category":"page"},{"location":"manual/installation/#Installation-of-the-package-PhyNEST","page":"Installation","title":"Installation of the package PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To use PhyNEST.jl, Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To install PhyNEST.jl run the following command in Julia prompt:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using Pkg\nPkg.add(\"PhyNEST\")","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] add PhyNEST","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Either command installs the package. To check the current version of installed PhyNEST.jl, run the following command:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] status PhyNEST","category":"page"},{"location":"manual/installation/#Update-PhyNEST","page":"Installation","title":"Update PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To update PhyNEST.jl to the most recent version, run:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Pkg.update(\"PhyNEST\")","category":"page"},{"location":"manual/installation/#Test-example","page":"Installation","title":"Test example","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Once the installation is complete, load PhyNEST.jl and run function greet() to see if the package has been added to the local machine. A greet message with current time should appear.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using PhyNEST\ngreet()","category":"page"},{"location":"manual/networkest/#Network-estimation","page":"Network estimation","title":"Network estimation","text":"","category":"section"},{"location":"manual/networkest/#Getting-ready-for-the-network-estimation","page":"Network estimation","title":"Getting ready for the network estimation","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"To start the network analysis, (1) sequence alignment file is required and (2) starting topology and (3) an outgroup must be specified.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"We use the PHYLIP alignment file in /example folder of PhyNEST package. The sequence alignment is parsed using function readPhylipFile! as shown here. ","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Then, we specify the starting tree (or network) stored in the HybridNetwork object using function readTopology() from the package PhyloNetworks. The starting tree can be either estimated from the data or randomly generated. Here, we use an arbitrarily generated tree topology.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"using PhyNEST\ndata = readPhylipFile!(\"n5h1_5k.txt\")\nstartingtree = readTopology(\"(5,(4,(3,(2,1))));\")","category":"page"},{"location":"manual/networkest/#Simulated-Annealing","page":"Network estimation","title":"Simulated Annealing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Function PhyNE!() executes the network analysis. By default, it uses simulated annealing algorithm to search the network space. PhyNE!() has a number of optional arguments (see here), but specifying the starting tree, data, and an outgroup and initiate the search as shown below: ","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netSA = PhyNE!(startingtree,data,\"5\")","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"To visualize the progress, the optional argument display=true can be used, which should print the following prompt, for example, at the initiation of the search.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST: Phylogenetic Network Estimation using SiTe patterns\nAnalysis start: 2022-11-07 at 21:21:22\nInput filename: n5h1_5k.txt\nNumber of sequences: 5 \nSequence length: 2500000\nStarting Topology: (5,(4,(3,(2,1))));\nOutgroup specified for rooting: 5\nNumber of maximum reticulation(s): 1\nThe maximum number of iterations for each optimization: 1000\nSearch algorithm selected: Simulated Annealing\nThe maximum number of steps during search: 100000\nAlpha: 0.8; Cons: 0.9\n\nInitiating 5 iterations...","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST.sa.log records the log throughout the searches. For each iteration (or independent run) using simulated annealing search, PhyNEST records k best networks searched along with other relevant information as shown below:","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"(1/5) Searching for the best network using the simulated annealing algorithm...\nStarting topology modified to (5,(4,(3,(2,1))));\nRunning 25 runs of burn-in...Complete \n(Cooling schedule)U=299136.023507182\n(Cooling schedule)beta=0.8999957157391261\nRank   Composite Likelihood    Network\n1\t2.88360994739e6         (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n2\t2.88360994739e6         (5:10.143246111976964,(4:5.128207274589672,((3:1.7832993566222581,(2)#H6:::0.5947475713189532):1.3303497322531534,(#H6:::0.40525242868104683,1:1.1961994949262897):1.9174495939491218):2.01455818571426):5.0150388373872925);\n3\t2.88362151446e6         (5:10.694491971616939,(4:5.43212278921928,((3:2.485548434584157,(#H6:::0.32538307301269104,2:6.112907457836508e-13):2.4855484345835457):1.2278094075673045,(1)#H6:::0.674616926987309):1.7187649470678181):5.262369182397659);\n4\t2.88362151446e6         (5:10.694492301887275,(4:5.43212302522998,((3:2.485548536348796,(2:3.738582439844855e-12,(1)#H6:::0.32538305925513744):2.4855485363450573):1.2278094281091798,#H6:::0.6746169407448626):1.7187650607720037):5.262369276657295);\n5\t2.88362693623e6         (5:10.725629067294516,(4:5.4494219219870095,(#H6:::0.6715226076537066,(1:2.508818645237626,(2:6.922689088873185e-13,(3)#H6:::0.32847739234629336):2.508818645236934):1.2147065868268707):1.7258966899225126):5.276207145307507);\n6\t2.88362693623e6         (5:10.725629357953252,(4:5.449422058791485,((1:2.508818707391581,(#H6:::0.32847740524951685,2:2.954945769153218e-13):2.5088187073912858):1.2147066492672614,(3)#H6:::0.6715225947504831):1.7258967021326423):5.276207299161768);\n7\t2.8837567506e6         (5:12.935996175051075,(#H6:::0.21325293589105276,((1:3.1629776123192412,(2:1.933312084245944,(3)#H6:::0.7867470641089472):1.2296655280732973):3.365071192479882,4:6.528048804799123):1.170591428480983):5.237355941770969);\n8\t2.88378761042e6         (5:13.131954265122538,((3:3.2109280989424116,(2:1.8718749279421563,(1)#H6:::0.7233614942899573):1.3390531710002553):3.5626956885379393,(4:6.773623787480175,#H6:::0.27663850571004267):1.758593271006248e-13):6.358330477642187);\n9\t2.88380850997e6         (5:13.030209792938253,((1:3.198422341899994,(2:1.8541586155663632,(3)#H6:::0.7273145740149164):1.3442637263336308):3.5191039964160598,(#H6:::0.2726854259850836,4:6.717526338316051):2.6645352591003757e-15):6.312683454622199);\n10\t2.88380850997e6         (5:13.030209981281956,(#H6:::0.27268542338357454,(4:6.71752636287232,(1:3.198422435121352,(2:1.8541587632941887,(3)#H6:::0.7273145766164255):1.3442636718271634):3.519103927750968):0.0):6.312683618409636);\nSpeciation times for some newicks may not have updated if estimates are weird (e.g., NaN).\nThe search terminated at step 360 and at 50th consecutive failures and (Cooling schedule)ci=922.9788432411981.\nSummary of each move:\nInsertion of reticulation edge: 1 proposed, 1 accepted.\nTail move of reticulation edge: 99 proposed, 37 accepted. \nHead move of reticulation edge: 101 proposed, 11 accepted.\nChange the direction of reticulation edge: 82 proposed, 19 accepted.\nDeletion of reticulation edge: 0 proposed, 0 accepted.\nNearest-neighbor interchange (NNI): 77 proposed, 9 accepted.\nOn the current topology, 60 moves were made, including 10 unsuccessful moves.\nTerminated because it reached the maximum number of failures (current nfail=50).\nThe best network found in this run: (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n-Log Composite Likelihood: 2.8836099473859877e6","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST.sa.out reports the maximum composite likelihood network found during the multiple runs of searches (see below for visualization). ","category":"page"},{"location":"manual/networkest/#Hill-climbing","page":"Network estimation","title":"Hill climbing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Network search using hill climbing algorithm can be initiated using an option hillclimbing=true.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netHC = PhyNE!(startingtree,data,\"5\",hillclimbing=true)","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Similar to the simulation annealing, it will create two output files named PhyNe.hc.log and PhyNe.hc.out Unlike the simulated annealing that provided k best networks per independent run, hill climbing yields a single maximum composite likelihood network for each run. The best network found at each run is logged in the output file PhyNEST.hc.log.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"(1/5) Searching for the best network using the hill-climbing algorithm...\nStarting topology modified to (5,(4,(3,(2,1))));\nThe search terminated at step 85 and at 75th consecutive failures.\nSummary of each move:\nInsertion of reticulation edge: 1 proposed, 1 accepted.\nTail move of reticulation edge: 26 proposed, 1 accepted. \nHead move of reticulation edge: 28 proposed, 0 accepted.\nChange the direction of reticulation edge: 8 proposed, 1 accepted.\nDeletion of reticulation edge: 0 proposed, 0 accepted.\nNearest-neighbor interchange (NNI): 22 proposed, 0 accepted.\nOn the current topology, 96 moves were made, including 21 unsuccessful moves.\nTerminated because it reached the maximum number of failures (current nfail=75).\nThe best network found in this run: (5:10.694492301887275,(4:5.43212302522998,((3:2.485548536348796,(2:3.738582439844855e-12,(1)#H6:::0.32538305925513744):2.4855485363450573):1.2278094281091798,#H6:::0.6746169407448626):1.7187650607720037):5.262369276657295);\n-Log Composite Likelihood: 2.883621514463323e6","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST.hc.out reports the maximum composite likelihood network found during the multiple runs of searches (see below for visualization). ","category":"page"},{"location":"manual/networkest/#Network-visualization","page":"Network estimation","title":"Network visualization","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"At the end of the search either using simulated annealing or hill climbing algorithm, the output file with an extension .out is created in the working directory that should look something similar to below:","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"The best network found from 5 runs using the simulated annealing algorithm\nMCL network: \n(5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,#H6:::0.4052526191450872):1.9174483823877784,((2)#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n    Dendroscope: \n(5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,#H6):1.9174483823877784,((2)#H6,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n    -Log Composite Likelihood: 2.8836099473859877e6.\nend","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Two extended Newick strings (See Cardona et al.,(2008)) represent the identical network topology, but formatted for visualization using the Julia package PhyloPlots.jl (top Newick) or Dendroscope 3 (bottom Newick).","category":"page"},{"location":"manual/others/#Others","page":"Others","title":"Others","text":"","category":"section"},{"location":"manual/others/#HyDe","page":"Others","title":"HyDe","text":"","category":"section"},{"location":"manual/others/#The-D-statistics","page":"Others","title":"The D-statistics","text":"","category":"section"},{"location":"manual/others/","page":"Others","title":"Others","text":"Manual TBA","category":"page"},{"location":"manual/quartet/#Quartet","page":"Quartet","title":"Quartet","text":"","category":"section"},{"location":"manual/quartet/#Reading-a-network","page":"Quartet","title":"Reading a network","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"An extended Newick string can be read by the function readTopology from PhyloNetworks. See Cardona et al.,(2008) to see the first description of the extended Newick formatted network. Here we use a five tip network:","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"PhyNEST uses the same the extended Newick format as in PhyloNetworks. Branch lengths can be specified using colon (:) as in a regular Newick string. For reticulation nodes, relevant information are specified in the order of ':length:bootstrap:gamma'. PhyNEST does not require branch lengths to be specified, however, if gamma is specified, it will be set as 0.5 .","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"using PhyNEST\nnetwork = readTopology(\"(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));\")","category":"page"},{"location":"manual/quartet/#Extract-quartet(s)","page":"Quartet","title":"Extract quartet(s)","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"PhyNE can extract quartets extracted from a tree or a network using the function extractQuartets.  ","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"printQuartets(network)","category":"page"},{"location":"manual/quartet/#True-site-pattern-probabilities-for-a-quartet","page":"Quartet","title":"True site pattern probabilities for a quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"See Chifman and Kubatko (2015) for more information.","category":"page"},{"location":"manual/quartet/#Symmetric-quartet","page":"Quartet","title":"Symmetric quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"True site pattern probabilities for the fifteen site pattern for a symmetric quartet under the JC69 model can be computed using the function TrueSitePatternSymm. We need to specify at least four parameters, three values of node ages, tau, in the order of [MRCA of species 1 and 2, MRCA of species 3 and 4, root age] and population size parameter theta. In the following block, we computed the probabilities for tau_1=10, tau_2=20, tau_3=50, and theta=00025. The 15 patterns appear in the order of [AAAA,AAAB,AABA,AABB,AABC,ABAA,ABAB,ABAC,ABBA,BAAA,ABBC,CABC,BACA,BCAA,ABCD]. ","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"We assure that the computation is on the right track by multiplying the weight for each site pattern, and sum all of them together to get 1.0.","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"symqProb=TrueSitePatternSymm(1.0,2.0,5.0,0.0025)\nweights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]\nsum(symqProb.*weights)","category":"page"},{"location":"manual/quartet/#Asymmetric-quartet","page":"Quartet","title":"Asymmetric quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"PhyNEST also computes true site pattern probabilities for the asymmetric quartet under the JC69 model using the function TrueSitePatternAsymm. The asymmetric quartet that we consider here can be written as (4,(3,(2,1))); in a Newick format. Similar to the symmetric quartet case above, three speciation times should be specified in the order of [MRCA of (1,2), MRCA of (1,2,3), root age], plus theta. ","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"asymqProb=TrueSitePatternAsymm(1.0,2.0,5.0,0.0025) \nweights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]\nsum(symqProb.*weights)","category":"page"},{"location":"manual/quartet/#Simulate-true-site-pattern-frequencies","page":"Quartet","title":"Simulate true site pattern frequencies","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"Function simspcounts generates the site-pattern frequencies for the symmetric and asymmetric quartets, modeled as a multinomial random variables assuming the observed sites are independent from each other. The length of sequence is 1000000 bp by default. Three speciation times for the quartet must be specified by the user. Theta and alpha are set as 0.0025 and 4/3, respectively, by default.","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"simSPsym=simspcounts(0,1.0,2.0,5.0,0.0025,4/3,1000000)\nsimSPasym=simspcounts(3,1.0,2.0,5.0,0.0025,4/3,1000000)","category":"page"},{"location":"manual/quartet/#Method-of-moment-estimator-of-branch-lengths","page":"Quartet","title":"Method-of-moment estimator of branch lengths","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function momentEstimat(type::Integer,spcounts::Array,theta::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"momEstsym=momentEstimat(0,simSPsym,0.0025)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"momEstasym=momentEstimat(3,simSPasym,0.0025)","category":"page"},{"location":"manual/quartet/#Estimating-theta","page":"Quartet","title":"Estimating theta","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function startTheta(q::Array{Nquartets, 1},net::HybridNetwork; lbound=0.00001::Float64,factor=2.0::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"datapath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"n5h1_5k.txt\");\nphydata = readPhylipFile!(datapath, showProgress=false)\ntheta=startTheta(network,phydata)","category":"page"},{"location":"#PhyNEST.jl","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Documentation for a Julia package PhyNEST.jl: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"page"},{"location":"#References","page":"PhyNEST.jl","title":"References","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Currently available on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.","category":"page"},{"location":"#Getting-help","page":"PhyNEST.jl","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Manual","page":"PhyNEST.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Pages = [\n    \"manual/installation.md\",\n    \"manual/input.md\",\n    \"manual/networkest.md\",\n    ]","category":"page"},{"location":"#Functions","page":"PhyNEST.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    extractNQuartets\n    extractQuartets\n    printQuartets\n    TrueSitePatternSymm\n    TrueSitePatternAsymm\n    simspcounts\n    PhyNE!","category":"page"},{"location":"#PhyNEST.readPhylipFile!","page":"PhyNEST.jl","title":"PhyNEST.readPhylipFile!","text":"readPhylipFile!(inputfile::AbstractString)\nreadPhylipFile!(inputfile::AbstractString; writecsv=false::Bool, csvname=\"\"::AbstractString, showProgress=true::Bool)\n\nTakes in, read, and parse the input phylip file. File name must be specified, in string format. There are a number of  optional arguments (see below).\n\nInput\n\ninputfile Name of phylip file as a AbstractString\nwritecsv  A Boolean arguemtn to write site pattern frequencies in CSV file (Default=false)\ncsvname   A string that will be name of the .csv file (Default=sitePatternCounts_inputfile.csv)\nshowProgress  A boolean argument for visualizing the progress of alignment parsing process (Default=true)\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readCheckPoint","page":"PhyNEST.jl","title":"PhyNEST.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file that is created every time function readPhylipFile! is complete.\n\nInput\n\nckpfile Name of checkpoint file as a AbstractString\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.writeSitePatternCounts","page":"PhyNEST.jl","title":"PhyNEST.writeSitePatternCounts","text":"writeSitePatternCounts(phylip::Phylip,csvname::AbstractString)\n\nExports the site pattern frequencies parsed frin the sequence alignment in a .csv file. The output is stored in the worksing directory unless specified, and named as sitePatternCounts_csvname.csv. \n\nInput\n\nphylip A Phylip object that is obtained using the function readPhylipFile!()\ncsvname  A string that specifies the name of the .csv file\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.extractNQuartets","page":"PhyNEST.jl","title":"PhyNEST.extractNQuartets","text":"extractNQuartets(net::HybridNetwork, phylip::Phylip)\n\nExtracts quartet information from a topology. Also stores the observed site patter frequency information for each quartet.\n\nInput\n\nnet   HybridNetwork object obtained using the function readTopology\nphylip    Phylip object obtained using the function readPhylipFile!\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.extractQuartets","page":"PhyNEST.jl","title":"PhyNEST.extractQuartets","text":"extractQuartets(net::HybridNetwork)\n\nExtracts quartet information from a topology.\n\nInput\n\nnet   HybridNetwork object obtained using the function readTopology\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.printQuartets","page":"PhyNEST.jl","title":"PhyNEST.printQuartets","text":"printQuartets(net::HybridNetwork)\n\nVisualizes quartet information extracted from a topology in the form of DataFrame.\n\nInput\n\nnet HybridNetwork object obtained using the function readTopology\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.TrueSitePatternSymm","page":"PhyNEST.jl","title":"PhyNEST.TrueSitePatternSymm","text":"TrueSitePatternSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nTrueSitePatternSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the symmetric quartet tree, ((1,2),(3,4));.  Three speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default.  See manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nInput\n\nmyt1 Speciation time for the common ancestor of species 1 and 2\nmyt2 Speciation time for the common ancestor of species 3 and 4\nmyt3 Root age\ntheta Effective population size parameter\nalpha 4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.TrueSitePatternAsymm","page":"PhyNEST.jl","title":"PhyNEST.TrueSitePatternAsymm","text":"TrueSitePatternAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nTrueSitePatternAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the asymmetric quartet tree, (1,(2,(3,4)));.  Three speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default.  See manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nInput\n\nmyt1 Speciation time for the common ancestor of species 3 and 4\nmyt2 Speciation time for the common ancestor of species 2, 3 and 4\nmyt3 Root age\ntheta Effective population size parameter\nalpha 4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.simspcounts","page":"PhyNEST.jl","title":"PhyNEST.simspcounts","text":"simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64)\nsimspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,length::Integer)\n\nThis function simulates the 15 site pattern probabilities for five quartet tree topologies,  one symmetric and four asymmetric quartets, which we categorize them into five types: \n\nType 0: ((1,2),(3,4));\nType 1: (1,((2,3),4));\nType 2: ((1,(2,3)),4);\nType 3: (1,(2,(3,4)));\nType 4: (((1,2),3),4);\n\nInput\n\ntype   Specify the type of a quartet (use integer)\nmyt1   Speciation time for the most recent internal tree node. Common ancestor of 1 and 2 in the symmetric case\nmyt2   Speciation time for the internal tree node closer to the root\nmyt3   Root age\ntheta  Effective population size parameter (default=0.01)\nalpha  Alpha (default=4/3)\nlength Sequence lengths (default=1000000)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.PhyNE!","page":"PhyNEST.jl","title":"PhyNEST.PhyNE!","text":"PhyNE!(startT::HybridNetwork,inputFile::Phylip,outgroup::String)\n\nEstimates the network or tree to fit observed site pattern frequencies stored in a Phylip object, using composite likelihood. A level-1 network is assumed. The search begins from topolgoy startT, which can be a tree or a network (with less than hmax reticulation nodes). Usually, startT is expected to be a tree estiamted from data, but can also be randomly generated. The topology is rooted using the outgroup provided. This must be identical to the outgroup sequence provided in the  sequence alignment. By default, PhyNE! will search the network space using simulated annealing assuming  hmax=1.\n\nThere are many optional arguments (values in parenthesis are defaults):\n\nhillclimbing(=false) Select hill climbing search\nhmax(=1) Maximum number of reticulations in the final network\nnruns(=10) Number of independent runs\nnniStartT=false \nmaxcount=1000 Number of steps to terminate the search\nnfail=75 Number of consecutive failures to terminate the search\n\nSimulated annealing stuff:\n\nburninn=25\nk=10\nProbabilityOfNotSelectingNode=9.5e-45\ncons=0.9\nalph=0.8\n\nOptimization stuff:\n\nNumIter=1000\nparamprint=false\n\nMiscellaneous:\n\ntimestamp=false\nfilename=\"\"\ndisplay=falses\n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNEST.jl","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    Phylip\n    nquartets\n    Nquartets","category":"page"},{"location":"#PhyNEST.Phylip","page":"PhyNEST.jl","title":"PhyNEST.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type. Sequence alignment with the following attributes:\n\nfilename    The actual filename of the input sequence alignment\n\ntime    Time it took to intake the filename in seconds\n\nnumtaxa The number of taxa appears on the first line of phylip\n\nseqleng The sequence length appears on the first line of phylip\n\nnametaxa    Names of each individual in the order of appearance in phylip\n\ncounttaxa   Numbering each individual in the order of appearance in phylip\n\nallquartet  All combinations of quartets using nametaxa\n\nindex   Just convert allquartet elements into arbitrary index numbers\n\nspcounts    Arrays of 15 site pattern frequencies \n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.nquartets","page":"PhyNEST.jl","title":"PhyNEST.nquartets","text":"nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a quartet extracted from a topology with the following attributes:\n\nnumber  List of individuals quartets in the order of i,j,k,l\n\ntquartet    List of quartets in the order of i,j,k,l\n\nindexNET    Index number for each quartet in the order appears in quartnet. \n\nmspcountsNET    We can directly use this counts for likelihood calculation.\n\nmrca    List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\n\nntau    Unique number of the taus in the tau used in branchlengths\n\nbranchlength    Branch length from the mrca to tau=0 in the order of i and j (ij),ik,il,jk,jl,and kl\n\nmomestlength    Branch length identified for each quartet given the mspcounts using moment estimator - Will get filled later because theta is required\n\nmombl   Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\n\nsymtype Type of each quartet. It can be either symmetric or asymmetric. Asymmetric (type 0) can have four possible topologies and we define type 1 as (i,((j,k),l)); type 2 as ((i,(j,k)),l); type 3 as (i,(j,(k,l))); and type 4 as (((i,j),k),l).\n\nlogLik  It's just a likelihood for the quartet\n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.Nquartets","page":"PhyNEST.jl","title":"PhyNEST.Nquartets","text":"Nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a tree extracted from a topology with the following attributes:\n\nleafname    Name of the leaves as appears on the newick\n\nleafnumber  Assigned number for each leaf when reading in the newick\n\ngamma   Gamma as appears on the newick, although irrelevant for our computation\n\nnquartet    An array of all quartets in the given topology\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNEST.jl","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    extractNQuartets\n    printQuartets\n    TrueSitePatternSymm\n    TrueSitePatternAsymm\n    PhyNE!\n    Phylip\n    nquartets\n    Nquartets","category":"page"}]
}
