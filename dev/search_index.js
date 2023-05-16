var documenterSearchIndex = {"docs":
[{"location":"manual/input/#Input","page":"Input","title":"Input","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Current version of PhyNEST.jl takes a relaxed, sequential, PHYLIP formatted DNA alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"A sample DNA alignment file sample_n5h1.txt is provided in the /example folder of the package PhyNEST.jl (or can be downloaded here). The alignment contains five sequences 1, 2, 3, 4, and 5 and is 1,000,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): (5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/input/#Parsing-DNA-alignment-data","page":"Input","title":"Parsing DNA alignment data","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Load PhyNEST.jl in Julia and move to the working directory that contains the alignment. The alignment is parsed with the function readPhylip(). The progress can be visualized with showProgress option that takes a boolean argument (set as false in this example for brevity). In the following command, we also set checkpoint option as true to create the .ckp file. Once the parsing the alignment is complete, a checkpoint file n5h1_5k.txt.ckp should appear in the directory. Note the .ckp file will have the same name as the alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"using PhyNEST\ndatapath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"sample_n5h1.txt\");\nphydata = readPhylip(datapath, showProgress=false, checkpoint=true)","category":"page"},{"location":"manual/input/#Reading-.ckp-file","page":"Input","title":"Reading .ckp file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"In case of multiple network analyses using the same alignment, a user can bypass the data parsing process by importing the .ckp file that is created after parsing the alignment for the first time. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Reading a checkpoint file using the function readCheckPoint() as shown below.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"ckppath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"sample_n5h1.txt.ckp\");\nckpdata = readCheckPoint(ckppath)","category":"page"},{"location":"manual/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual/installation/#Julia-installation","page":"Installation","title":"Julia installation","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single programming language. Current release of Julia is available here. Julia can be credited by citing the following paper:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65–98. doi: 10.1137/141000671.","category":"page"},{"location":"manual/installation/#Installation-of-the-package-PhyNEST","page":"Installation","title":"Installation of the package PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To use PhyNEST.jl, Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To install PhyNEST.jl run the following command in Julia prompt:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using Pkg\nPkg.add(\"PhyNEST\")","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] add PhyNEST","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Either command installs the package. To check the current version of installed PhyNEST.jl, run the following command:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] status PhyNEST","category":"page"},{"location":"manual/installation/#Update-PhyNEST","page":"Installation","title":"Update PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To update PhyNEST.jl to the most recent version, run:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Pkg.update(\"PhyNEST\")","category":"page"},{"location":"manual/installation/#Test-example","page":"Installation","title":"Test example","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Once the installation is complete, load PhyNEST.jl and run function greet() to see if the package has been added to the local machine. A greet message with current time should appear.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using PhyNEST\ngreet()","category":"page"},{"location":"manual/networkest/#Network-estimation","page":"Network estimation","title":"Network estimation","text":"","category":"section"},{"location":"manual/networkest/#Getting-ready-for-the-network-estimation","page":"Network estimation","title":"Getting ready for the network estimation","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"To start the network analysis, we need (1) starting topology, (2) sequence data and (3) an outgroup taxa.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"We use the PHYLIP alignment file in /example folder of PhyNEST package. The sequence alignment is parsed using function readPhyli as shown here. ","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Then, we specify the starting tree (or network) stored in the HybridNetwork object using function readTopology() from the package PhyloNetworks. The starting tree can be either estimated from the data or randomly generated. Here, we use an arbitrary tree topology.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"using PhyNEST\ndata = readPhylipFile!(\"n5h1_5k.txt\")\nstartingtree = readTopology(\"(5,(4,(3,(2,1))));\")","category":"page"},{"location":"manual/networkest/#Hill-climbing","page":"Network estimation","title":"Hill climbing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Network search using hill climbing algorithm is set as default","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netHC = phyne!(startingtree,data,\"5\")","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"The following will be displayed STDOUT.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST: Phylogenetic Network Estimation using SiTe patterns\nAnalysis start: 2022-11-07 at 21:21:22\nInput filename: n5h1_5k.txt\nNumber of sequences: 5 \nSequence length: 2500000\nStarting Topology: (5,(4,(3,(2,1))));\nOutgroup specified for rooting: 5\nNumber of maximum reticulation(s): 1\nThe maximum number of iterations for each optimization: 1000\nSearch algorithm selected: Simulated Annealing\nThe maximum number of steps during search: 100000\nAlpha: 0.8; Cons: 0.9\n\nInitiating 5 iterations...","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"The analysis will create an output file named PhyNEST_hc.log. At each independent run, the analysis yields a single maximum composite likelihood network found.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"(1/15) Searching for the best network using the hill climbing algorithm...The search terminated at step 180 and at 100th consecutive failures.\nSummary of each move:\nInsertion of reticulation edge: 2 proposed, 2 accepted.\nTail move of reticulation edge: 56 proposed, 2 accepted. \nHead move of reticulation edge: 52 proposed, 0 accepted.\nChange the direction of reticulation edge: 27 proposed, 1 accepted.\nDeletion of reticulation edge: 0 proposed, 0 accepted.\nNearest-neighbor interchange (NNI): 43 proposed, 1 accepted.\nOn the current topology, 115 moves were made, including 15 unsuccessful moves.\nTerminated because it reached the maximum number of failures (current maximum_number_of_failures=100).\n(1/15) Estimated topology in this run: (7:5.858901865538142,(((((1:0.4895402649806814,(2:1.263190463023415)#H12:0.0::0.4954631266401886):0.7736501980427335,(3:0.4872669126554865,#H12:0.0::0.5045368733598115):0.7759235503679285):1.2442147399714147,4:2.5074052029948297):0.0,(5:0.0)#H9:0.0::0.4084401391767469):0.0,(6:0.325378144365268,#H9:0.0::0.5915598608232531):2.182027058629562):3.3514966625433122);\n(1/15) Composite likelihood of the estimated topology in this run: 6.405284917395862e7","category":"page"},{"location":"manual/networkest/#Simulated-annealing","page":"Network estimation","title":"Simulated annealing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Function phyne!() can search the network space using simulated annealing algorithm by setting `dohillclimbing=false'.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netSA = phyne!(startingtree,data,\"5\",do_hill_climbing=false)","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNEST_sa.log records the log throughout the searches. For each iteration (or independent run) using simulated annealing search, PhyNEST records k best networks searched along with other relevant information as shown below:","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"(1/5) Searching for the best network using the simulated annealing algorithm...\nStarting topology modified to (5,(4,(3,(2,1))));\nRunning 25 runs of burn-in...Complete \n(Cooling schedule)U=299136.023507182\n(Cooling schedule)beta=0.8999957157391261\nRank   Composite Likelihood    Network\n1\t2.88360994739e6         (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n2\t2.88360994739e6         (5:10.143246111976964,(4:5.128207274589672,((3:1.7832993566222581,(2)#H6:::0.5947475713189532):1.3303497322531534,(#H6:::0.40525242868104683,1:1.1961994949262897):1.9174495939491218):2.01455818571426):5.0150388373872925);\n3\t2.88362151446e6         (5:10.694491971616939,(4:5.43212278921928,((3:2.485548434584157,(#H6:::0.32538307301269104,2:6.112907457836508e-13):2.4855484345835457):1.2278094075673045,(1)#H6:::0.674616926987309):1.7187649470678181):5.262369182397659);\n4\t2.88362151446e6         (5:10.694492301887275,(4:5.43212302522998,((3:2.485548536348796,(2:3.738582439844855e-12,(1)#H6:::0.32538305925513744):2.4855485363450573):1.2278094281091798,#H6:::0.6746169407448626):1.7187650607720037):5.262369276657295);\n5\t2.88362693623e6         (5:10.725629067294516,(4:5.4494219219870095,(#H6:::0.6715226076537066,(1:2.508818645237626,(2:6.922689088873185e-13,(3)#H6:::0.32847739234629336):2.508818645236934):1.2147065868268707):1.7258966899225126):5.276207145307507);\n6\t2.88362693623e6         (5:10.725629357953252,(4:5.449422058791485,((1:2.508818707391581,(#H6:::0.32847740524951685,2:2.954945769153218e-13):2.5088187073912858):1.2147066492672614,(3)#H6:::0.6715225947504831):1.7258967021326423):5.276207299161768);\n7\t2.8837567506e6         (5:12.935996175051075,(#H6:::0.21325293589105276,((1:3.1629776123192412,(2:1.933312084245944,(3)#H6:::0.7867470641089472):1.2296655280732973):3.365071192479882,4:6.528048804799123):1.170591428480983):5.237355941770969);\n8\t2.88378761042e6         (5:13.131954265122538,((3:3.2109280989424116,(2:1.8718749279421563,(1)#H6:::0.7233614942899573):1.3390531710002553):3.5626956885379393,(4:6.773623787480175,#H6:::0.27663850571004267):1.758593271006248e-13):6.358330477642187);\n9\t2.88380850997e6         (5:13.030209792938253,((1:3.198422341899994,(2:1.8541586155663632,(3)#H6:::0.7273145740149164):1.3442637263336308):3.5191039964160598,(#H6:::0.2726854259850836,4:6.717526338316051):2.6645352591003757e-15):6.312683454622199);\n10\t2.88380850997e6         (5:13.030209981281956,(#H6:::0.27268542338357454,(4:6.71752636287232,(1:3.198422435121352,(2:1.8541587632941887,(3)#H6:::0.7273145766164255):1.3442636718271634):3.519103927750968):0.0):6.312683618409636);\nSpeciation times for some newicks may not have updated if estimates are weird (e.g., NaN).\nThe search terminated at step 360 and at 50th consecutive failures and (Cooling schedule)ci=922.9788432411981.\nSummary of each move:\nInsertion of reticulation edge: 1 proposed, 1 accepted.\nTail move of reticulation edge: 99 proposed, 37 accepted. \nHead move of reticulation edge: 101 proposed, 11 accepted.\nChange the direction of reticulation edge: 82 proposed, 19 accepted.\nDeletion of reticulation edge: 0 proposed, 0 accepted.\nNearest-neighbor interchange (NNI): 77 proposed, 9 accepted.\nOn the current topology, 60 moves were made, including 10 unsuccessful moves.\nTerminated because it reached the maximum number of failures (current nfail=50).\nThe best network found in this run: (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);\n-Log Composite Likelihood: 2.8836099473859877e6","category":"page"},{"location":"#PhyNEST.jl","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Documentation for a Julia package PhyNEST.jl: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"page"},{"location":"#References","page":"PhyNEST.jl","title":"References","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.","category":"page"},{"location":"#Getting-help","page":"PhyNEST.jl","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Manual","page":"PhyNEST.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Pages = [\n    \"manual/installation.md\",\n    \"manual/input.md\",\n    \"manual/networkest.md\",\n    ]","category":"page"},{"location":"#Functions","page":"PhyNEST.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    greet\n    readPhylip\n    show_sp\n    write_sp\n    storeCheckPoint\n    readCheckPoint\n    GetTrueProbsSymm\n    GetTrueProbsAsymm\n    GetTrueProbsAsymmTypes\n    sim_sp_counts\n    get_quartets\n    printQuartets\n    get_start_theta\n    do_optimization\n    hill_climbing\n    simulated_annealing\n    phyne!","category":"page"},{"location":"#PhyNEST.greet","page":"PhyNEST.jl","title":"PhyNEST.greet","text":"greet()\n\nDisplays a simple greet message with citation information. No input argument is needed. \n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readPhylip","page":"PhyNEST.jl","title":"PhyNEST.readPhylip","text":"readPhylip(inputfile::AbstractString)\n\nImport, read, and parse the input phylip file. File name must be specified as string. \nBy default, .csv and .ckp files are NOT produced. See below for the optional arguments.\\ \n\nTo save the observed quartet site pattern frequencies in a .csv file, set the optional argument\nwritecsv=true. The .csv file will be stored in the working directory, with a filename \n\"sitePatternCounts_inputfile.csv\". If a user would like to modify the filename, it can be done\nby providing preferred name using the optional argument csvname=\"preferred_file_name.csv\". \n\n\nInput\n\nMandatory\n\ninputfile     Name of the phylip file as a String [mandatory]\\ \n\nOptional\n\nwritecsv       (default=false) A Boolean arguement to write site pattern frequencies in CSV file\n\ncsvname        (default=sitePatternCounts_inputfile.csv) A string that will be name of the .csv file \n\nshowProgress   (default=true) A boolean argument for visualizing the progress of alignment parsing process\n\ntdigits        (default=3) Number of decimal points to record timme taken to parse the file in seconds\n\ncheckpoint     (default=false) A boolean argument to store Phylip object as a .ckp file. (Warning: .ckp can be large).\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.show_sp","page":"PhyNEST.jl","title":"PhyNEST.show_sp","text":"show_sp(p::Phylip)\n\nPretty name for displaying all quartet and the corresponding site pattern frequencies on screen in the DataFrame format. \nIt may result in excessively long table when there are many sequences in the input file.\n\nInput\n\np     Phylip object [mandatory]\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.write_sp","page":"PhyNEST.jl","title":"PhyNEST.write_sp","text":"write_sp(p::Phylip)\n\nPretty name for exporting all quartet and the corresponding site pattern frequencies in a .csv format. \n\n\nInput\n\np         Phylip object [mandatory] csvname   Filename for the .csv output can be given, otherwise we use PhyNEST_sp.csv by default.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.storeCheckPoint","page":"PhyNEST.jl","title":"PhyNEST.storeCheckPoint","text":"storeCheckPoint(p::Phylip)\n\nStores the phylip object into a .ckp file so a user does not have to repeat the phylip parsing again if working with the \nsame dataset next time. By default it creates a .ckp file that has the same filename as the input phylip file.\nFile size can get pretty large... \n\nInput\n\np         Phylip object [mandatory]\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readCheckPoint","page":"PhyNEST.jl","title":"PhyNEST.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file and creates a filled in phylip object.\n\nInput\n\nckpfile   Name of the checkpoint file\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsSymm","page":"PhyNEST.jl","title":"PhyNEST.GetTrueProbsSymm","text":"GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nGetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the symmetric quartet tree, ((1,2),(3,4));. \nThree speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not \nan essential argument for the function and if not provided, it isssumed to be 4/3 by default. \nSee manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\n\nInput\n\nmyt1      Speciation time for the common ancestor of species 1 and 2\nmyt2      Speciation time for the common ancestor of species 3 and 4\nmyt3      Root node age\ntheta     Effective population size parameter\nalpha     4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsAsymm","page":"PhyNEST.jl","title":"PhyNEST.GetTrueProbsAsymm","text":"GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nGetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the asymmetric quartet tree, (1,(2,(3,4)));. \nThree speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not\nan essential argument for the function and if not provided, it isssumed to be 4/3 by default. \nSee manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\n\nInput\n\nmyt1      Speciation time for the common ancestor of species 3 and 4\nmyt2      Speciation time for the common ancestor of species 2, 3 and 4\nmyt3      Root node age\ntheta     Effective population size parameter\nalpha     4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsAsymmTypes","page":"PhyNEST.jl","title":"PhyNEST.GetTrueProbsAsymmTypes","text":"GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nGetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for any of the four the asymmetric quartet trees: \n\n\nType 1: (i,((j,k),l));\nType 2: ((i,(j,k)),l); \nType 3: (i,(j,(k,l))); \nType 4: (((i,j),k),l).\n\n\nThree speciation times (or node ages) in coalescent unit and theta must be provided. \nAlhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default. \nSee manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\n\nInput\n\nmyt1      Speciation time for the common ancestor of species 3 and 4\nmyt2      Speciation time for the common ancestor of species 2, 3 and 4\nmyt3      Root node age\ntheta     Effective population size parameter\nalpha     4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sim_sp_counts","page":"PhyNEST.jl","title":"PhyNEST.sim_sp_counts","text":"sim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64)\nsim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,length::Integer)\n\nThis function simulates the 15 site pattern probabilities for five quartet tree topologies, \none symmetric and four asymmetric quartets: \n\nType 0: ((i,j),(k,l));\nType 1: (i,((j,k),l));\nType 2: ((i,(j,k)),l); \nType 3: (i,(j,(k,l))); \nType 4: (((i,j),k),l).\n\n\nInput\n\ntype   Specify the type of a quartet (use integer)\nmyt1   Speciation time for the most recent internal tree node. Common ancestor of 1 and 2 in the symmetric case\nmyt2   Speciation time for the internal tree node closer to the root\nmyt3   Root node age\ntheta  Effective population size parameter (default=0.01)\nalpha  Alpha (default=4/3)\nlength Sequence lengths (default=1000000)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.get_quartets","page":"PhyNEST.jl","title":"PhyNEST.get_quartets","text":"get_quartets(net::HybridNetwork; round_gamma_digits=5::Int64)\nget_quartets(net::HybridNetwork, p::Phylip; round_gamma_digits=5::Int64, default_theta=0.0001::Float64)\n\nExtracts information for every quartet extractable from the topology and store it as a Network object. \nWhen the Phylip object is also provided, we fill more slots including the theta, branch lengths for \neac quartet using the moment estimator, average branch lengths using the estimated branch lengths.\n\nInput\n\nnet                   PhyloNetwork.HybridNetwork object for a topology [mandatory]\n\np                     Phylip object\n\nround_gamma_digits    Number of digits to show for the inhertiance probability in each quartet\n\n`default_theta          Theta is being estimated when both topology and data are given, however, \n                      when the estimated theta is beyond the expected range (i.e., 0.00001 < theta < 0.1)\n                      we use user-specified theta (default=0.0001) to estimate branch lengths using \n                      the moment estimator.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.get_start_theta","page":"PhyNEST.jl","title":"PhyNEST.get_start_theta","text":"get_start_theta(N::Network; \n                lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)\n\nUsing the Network object that contains information on quartets extracted from the topology and\nsite pattern frequencies from the data, estimates the `reasonable' theta value that can be\nused as the starting point for the composite likelihood optimization. \nApproximate interval where the true theta may lie is estimated through two procedures:\nsomewhat loose interval that keeps the branch lengths positive from the moment estimator,\nand then further tightens the interval using the golden section seach.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.do_optimization","page":"PhyNEST.jl","title":"PhyNEST.do_optimization","text":"do_optimization(net::HybridNetwork, p::Phylip;\n                    lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64,\n                    number_of_itera=1000::Int64)\n\nOptimizes the composite likelihood of the network topology given the data.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.hill_climbing","page":"PhyNEST.jl","title":"PhyNEST.hill_climbing","text":"hill_climbing(starting_topology::HybridNetwork,p::Phylip,outgroup::String;\n    hmax::Integer,\n    maximum_number_of_steps::Integer,\n    maximum_number_of_failures::Integer,\n    number_of_itera::Integer)\n\nExecutes hill climbing algorithm for the network space search.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.simulated_annealing","page":"PhyNEST.jl","title":"PhyNEST.simulated_annealing","text":"simulated_annealing(starting_topology::HybridNetwork, p::Phylip, outgroup::String, \n    hmax::Integer, \n    number_of_burn_in::Int64,\n    maximum_number_of_steps::Integer,\n    maximum_number_of_failures::Integer,\n    number_of_itera::Integer,\n    k::Integer,\n    cons::Float64,\n    alph::Float64\n    )\n\nExecutes simulated algorithm for the network space search.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.phyne!","page":"PhyNEST.jl","title":"PhyNEST.phyne!","text":"phyne!(starting_topology::HybridNetwork,p::Phylip,outgroup::String;\n    hmax=1::Integer,\n    maximum_number_of_steps=1000::Integer,\n    maximum_number_of_failures=100::Integer,\n    number_of_itera=1000::Integer,\n    number_of_runs=5::Integer,\n    do_hill_climbing=true::Bool,\n    number_of_burn_in=25::Integer,\n    k=10::Integer,\n    cons=0.9::Float64,\n    alph=0.8::Float64\n    )\n\nphyne! is fine.\n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNEST.jl","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    Phylip\n    Network\n    quartets","category":"page"},{"location":"#PhyNEST.Phylip","page":"PhyNEST.jl","title":"PhyNEST.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type with the following attributes:\n\nfilename      Name of the input phylip alignment file\n\ntime          Time taken for parsing the input in seconds\n\nnumtaxa       Number of taxa given at the first line of the input file\n\nseqleng       Sequence length given at the first line of the input file\n\nnametaxa      Sequence names given in the input file\n\ncounttaxa     Unique integer identifier given to each individual in the order of appearance in the input file\n\nallquartet    All combinations of quartets for counttaxa\n\nindex         Just convert allquartet elements into arbitrary index numbers\n\nspcounts      Arrays of 15 site pattern frequencies for each quartet in allquartet\n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.Network","page":"PhyNEST.jl","title":"PhyNEST.Network","text":"Network\n\nSubtype of abstract Quartet type with the following attributes:\n\nleafname    Name of the leaves in the order of appearance in the input newick\n\nleafnumber  Assigned number for each leaf when reading in the input newick\n\ngamma       Gamma as appears on the newick, although irrelevant for our computation\n\ntheta       Estimated theta using the site pattern frequencies. If Network object is\n              created without Phylip, theta is not estimated byt default theta=0.001 is displayed.\n\nquartet     An array of all quartets in the given topology\n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.quartets","page":"PhyNEST.jl","title":"PhyNEST.quartets","text":"quartets\n\nSubtype of abstract Quartet type with the following attributes:\n\nnumber        List of individuals quartets in the order of i,j,k,l\ndisplayed_treenth displayed tree that the quartet was extracted from\nquartet       List of quartets in the order of i,j,k,l using the leaf numbers in HybridNetwork\ntquartet      List of quartets in the order of i,j,k,l using the leaf numbers in Phylip\ngamma         Inhertiance probability information provided in HybridNetwork\nmspcountsNET  We can directly use this counts for likelihood calculation.\nmrca          List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\nntau          Unique number of the taus in the tau used in branchlengths\nmomestlength  Branch length identified for each quartet given the mspcounts using moment estimator\naverage_mom_est_bl         Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\nsymtype       Type of each quartet. It can be either symmetric (type 0) or asymmetric. Asymmetric quartets have four possible topologies:\n\nType 1: (i,((j,k),l));\nType 2: ((i,(j,k)),l); \nType 3: (i,(j,(k,l))); \nType 4: (((i,j),k),l).\n\n\nlogLik        It's just a negative likelihood for the quartet in interest\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNEST.jl","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    readPhylip\n    readPhylipFile\n    PhylipFileInfo\n    getUniqueQuartets\n    sitePatternCounts\n    spRearrange\n    show_sp\n    sitePatternsToDF\n    write_sp\n    writeSitePatternCounts\n    storeCheckPoint\n    readCheckPoint\n    get_quartets\n    get_leaf_name\n    get_leaf_number\n    list_all_quartets\n    GetChild\n    GetParent\n    childnode\n    parentnode\n    mrca\n    get_most_recent_common_ancestors\n    identify_sym_type\n    tauNum\n    backtauNum\n    get_unique_tau_labels\n    dictionary_phylip\n    dictionary_topology\n    dictionary_combined\n    get_matching_quartet\n    transfer_site_pattern_frequencies\n    get_gamma\n    printQuartets\n    GetTrueProbsSymm\n    GetTrueProbsAsymm\n    GetTrueProbsAsymmTypes\n    GetTrueProbsNetTypes\n    sim_sp_counts\n    do_optimization\n    get_start_gamma_info\n    get_parent_child\n    get_initial_taus\n    backTransform\n    backtransform_parameters\n    update_topology\n    hill_climbing\n    burn_in\n    simulated_annealing\n    initiate_search\n    phyne!\n    get_average\n    greet\n    checkDEBUG\n    Dstat\n    HyDe","category":"page"}]
}
