var documenterSearchIndex = {"docs":
[{"location":"manual/input/#Input","page":"Input","title":"Input","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Current version of PhyNE takes a relaxed, sequential, PHYLIP formatted DNA alignment file. Extensions can be either .phy or .txt. For more information on PHYLIP, see here. After the PHYLIP file has been parsed, PhyNE creates a checkpoint file with the extension .ckp, which can be called later to prevent repetition of parsing the same data and save time. Moreover, the observed site pattern frequencies of every quartet from the alignment can be exported into a .csv file, for a user to explore for other purposes (e.g., manually computation of Patterson's D-statistic etc.). The following sections explain above procedures in detail.","category":"page"},{"location":"manual/input/#Parsing-DNA-alignment-data","page":"Input","title":"Parsing DNA alignment data","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"A sample DNA alignment file called n5h1_5k.txt is given in the example folder of the package PhyNE (or can be downloaded here). This alignment is 750,000 base pairs (bp) long and contains five sequences, named as  [1, 2, 3, 4, 5]. The alignment is generated using ms and seq-gen. The true relationship of the five taxa is written in the extended Newick format as (branch lengths omittted): ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"We begin with loading the package PhyNE in Julia.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"using PhyNEST","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Assuming that you are working in the directory that contains the alignment file, the alignment file can be read using the function readPhylipFile!(inputfile). Since the variable inputfile is expected to be a string, the file name should be in the quotation marks. The path of the file can be set as shown below or simply the name of the alignment file can be specified (e.g., readPhylipFile!(\"n5h1_5k.txt\")). A checkpoint file with an extension .ckp is automatically created in the working directory upon completion of parsing the alignment. When the boolean option showProgress is set as true, PhyNE will show the progress bar during the data parsing process. Here, we set it as false for brevity.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"datapath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"n5h1_5k.txt\");\nphydata = readPhylipFile!(datapath, showProgress=false)","category":"page"},{"location":"manual/input/#Reading-.ckp-file","page":"Input","title":"Reading .ckp file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Every time a sequence alignment is parsed, PhyNE creates a checkpoint file with the extension .ckp. Note the .ckp file will have the same name as the alignment file. In case multiple network analyses using the same data is planned, a user can bypass the (potentially time-consuming) data parsing process by using the function readCheckPoint(ckpfile) instead of the function readPhylipFile!(inputfile). ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"ckppath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"n5h1_5k.txt.ckp\");\nckpdata = readCheckPoint(ckppath)","category":"page"},{"location":"manual/input/#Exporting-as-.csv-file","page":"Input","title":"Exporting as .csv file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"The observed quartet site pattern frequencies from the data can be exported as a .csv formatted file. This can be done when parsing the alignment using the function readPhylipFile! by setting the option writecsv=true. It will create a .csv file in the working directory with the filename sitePatternCounts_$inputfile.csv. If a user would like to save the .csv file in a different name, the file name can be specified in the option csvname. In the below example, we set the filename as \"tutorial_n5h1\" for the .csv file. The .csv file can be visualized using function readCSVFile(csvfile). ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"If a user forgot to set writecsv=true, the function writeSitePatternCounts can be used to export. See here for more information.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"phydata = readPhylipFile!(datapath, showProgress=false,\n                        writecsv=true, csvname=\"tutorial_n5h1\")\ncsvpath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"tutorial_n5h1.csv\");\ncsvdf = readCSVFile(csvpath)","category":"page"},{"location":"manual/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual/installation/#Julia-installation","page":"Installation","title":"Julia installation","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single language. Current release of Julia is available here. Julia can be creditted by citing the following paper:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Julia: A Fresh Approach to Numerical Computing. Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) SIAM Review, 59: 65–98. doi: 10.1137/141000671.","category":"page"},{"location":"manual/installation/#Installation-of-the-package-PhyNE","page":"Installation","title":"Installation of the package PhyNE","text":"","category":"section"},{"location":"manual/installation/#Test-example","page":"Installation","title":"Test example","text":"","category":"section"},{"location":"manual/networkest/#Network-estimation","page":"Network estimation","title":"Network estimation","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"After parsing the input sequence alignment into quartet site pattern frequency data, we can estimate the network using a starting tree/network, that can be either randomly generated or estimated from the data. ","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"startingtree = readTopology(\"(5,(4,(3,(2,1))));\")\ndata=readPhylipFile!(\"phylipfile.phy\")","category":"page"},{"location":"manual/networkest/#Hill-climbing","page":"Network estimation","title":"Hill climbing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"Network search using hill climbing algorithm can be initiated using:","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netHC = PhyNE!(startingtree,data,\"outgroup\",hillclimbing=true)","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"To visualize the progress, an option display=true can be used, and the following prompt should be printed on the screen that summarizes some information of the data and the searching condition.","category":"page"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"PhyNe: Estimating Maximum Pseudolikelihood Phylogenetic Network\n╔═══╦╗─────╔═╗─╔╗\n║╔═╗║║─────║║╚╗║║\n║╚═╝║╚═╦╗─╔╣╔╗╚╝╠══╗\n║╔══╣╔╗║║─║║║╚╗║║║═╣\n║║──║║║║╚═╝║║─║║║║═╣\n╚╝──╚╝╚╩═╗╔╩╝─╚═╩══╝\n───────╔═╝║ \n───────╚══╝   \nAnalysis start: 2022-10-07 at 16:08:02\nInput filename: n5h1_5k.txt\nNumber of sequences: 5 \nSequence length: 2500000\nStarting Topology: (5,(4,(3,(2,1))));\nOutgroup specified for rooting: 1\nNumber of maximum reticulation(s): 1\nThe maximum number of iterations for each optimization: 1000\nSearch algorithm selected: Hill-climbing\nThe maximum number of steps during search: 50","category":"page"},{"location":"manual/networkest/#Simulated-Annealing","page":"Network estimation","title":"Simulated Annealing","text":"","category":"section"},{"location":"manual/networkest/","page":"Network estimation","title":"Network estimation","text":"netSA = PhyNE!(startingtree,data,\"outgroup\")","category":"page"},{"location":"manual/others/#Others","page":"Others","title":"Others","text":"","category":"section"},{"location":"manual/others/#HyDe","page":"Others","title":"HyDe","text":"","category":"section"},{"location":"manual/others/#The-D-statistics","page":"Others","title":"The D-statistics","text":"","category":"section"},{"location":"manual/others/","page":"Others","title":"Others","text":"Manual TBA","category":"page"},{"location":"manual/quartet/#Quartet","page":"Quartet","title":"Quartet","text":"","category":"section"},{"location":"manual/quartet/#Reading-a-network","page":"Quartet","title":"Reading a network","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"An extended Newick string can be read by the function readTopology from PhyloNetworks. See Cardona et al.,(2008) to see the first description of the extended Newick formatted network. Here we use a five tip network:","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"PhyNE uses the same the extended Newick format as in PhyloNetworks. Branch lengths can be specified using colon (:) as in a regular Newick string. For reticulation nodes, relevant information are specified in the order of ':length:bootstrap:gamma'. PhyNE does not require branch lengths to be specified, however, if gamma is specified, it will be set as 0.5 .","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"using PhyNEST\nnetwork = readTopology(\"(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));\")","category":"page"},{"location":"manual/quartet/#Extract-quartet(s)","page":"Quartet","title":"Extract quartet(s)","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"PhyNE can extract quartets extracted from a tree or a network using the function extractQuartets.  ","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"printQuartets(network)","category":"page"},{"location":"manual/quartet/#True-site-pattern-probabilities-for-a-quartet","page":"Quartet","title":"True site pattern probabilities for a quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"See Chifman and Kubatko (2015) for more information.","category":"page"},{"location":"manual/quartet/#Symmetric-quartet","page":"Quartet","title":"Symmetric quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"True probabilities for the fifteen site pattern for a symmetric quartet can be computed using the function GetTrueProbsSymm. We need to specify five parameters, tau","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"symqProb=GetTrueProbsSymm(1.0,2.0,5.0,0.0025,4/3)","category":"page"},{"location":"manual/quartet/#Asymmetric-quartet","page":"Quartet","title":"Asymmetric quartet","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"asymqProb=GetTrueProbsAsymm(1.0,2.0,5.0,0.0025,4/3) ","category":"page"},{"location":"manual/quartet/#Simulate-true-site-pattern-frequencies","page":"Quartet","title":"Simulate true site pattern frequencies","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"simSPsym=simspcounts(0,1.0,2.0,5.0,0.0025,4/3,1000000)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"simSPasym=simspcounts(3,1.0,2.0,5.0,0.0025,4/3,1000000)","category":"page"},{"location":"manual/quartet/#Method-of-moment-estimator-of-branch-lengths","page":"Quartet","title":"Method-of-moment estimator of branch lengths","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function momentEstimat(type::Integer,spcounts::Array,theta::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"momEstsym=momentEstimat(0,simSPsym,0.0025)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"momEstasym=momentEstimat(3,simSPasym,0.0025)","category":"page"},{"location":"manual/quartet/#Estimating-theta","page":"Quartet","title":"Estimating theta","text":"","category":"section"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"function startTheta(q::Array{Nquartets, 1},net::HybridNetwork; lbound=0.00001::Float64,factor=2.0::Float64)","category":"page"},{"location":"manual/quartet/","page":"Quartet","title":"Quartet","text":"datapath = joinpath(dirname(pathof(PhyNE)), \"..\",\"example\",\"n5h1_5k.txt\");\nphydata = readPhylipFile!(datapath, showProgress=false)\ntheta=startTheta(network,phydata)","category":"page"},{"location":"#PhyNEST.jl","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Documentation for a julia package PhyNEST.jl: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"page"},{"location":"#References","page":"PhyNEST.jl","title":"References","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"TBA","category":"page"},{"location":"#Getting-help","page":"PhyNEST.jl","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Manual","page":"PhyNEST.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Pages = [\n    \"manual/installation.md\",\n    \"manual/input.md\",\n    \"manual/networkest.md\",\n    ]","category":"page"},{"location":"#Functions","page":"PhyNEST.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    extractNQuartets\n    extractQuartets\n    printQuartets\n    TrueSitePatternSymm\n    TrueSitePatternAsymm\n    PhyNE!","category":"page"},{"location":"#PhyNEST.readPhylipFile!","page":"PhyNEST.jl","title":"PhyNEST.readPhylipFile!","text":"readPhylipFile!(inputfile::AbstractString)\nreadPhylipFile!(inputfile::AbstractString; writecsv=false::Bool, csvname=\"\"::AbstractString, showProgress=true::Bool)\n\nTakes in, read, and parse the input phylip file. File name must be specified, in string format. There are a number of  optional arguments (see below).\n\nInput\n\ninputfile Name of phylip file as a AbstractString\nwritecsv  A Boolean arguemtn to write site pattern frequencies in CSV file (Default=false)\ncsvname   A string that will be name of the .csv file (Default=sitePatternCounts_inputfile.csv)\nshowProgress  A boolean argument for visualizing the progress of alignment parsing process (Default=true)\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readCheckPoint","page":"PhyNEST.jl","title":"PhyNEST.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file that is created every time function readPhylipFile! is complete.\n\nInput\n\nckpfile Name of checkpoint file as a AbstractString\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.writeSitePatternCounts","page":"PhyNEST.jl","title":"PhyNEST.writeSitePatternCounts","text":"writeSitePatternCounts(phylip::Phylip,csvname::AbstractString)\n\nExports the site pattern frequencies parsed frin the sequence alignment in a .csv file. The output is stored in the worksing directory unless specified, and named as sitePatternCounts_csvname.csv. \n\nInput\n\nphylip A Phylip object that is obtained using the function readPhylipFile!()\ncsvname  A string that specifies the name of the .csv file\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.extractNQuartets","page":"PhyNEST.jl","title":"PhyNEST.extractNQuartets","text":"extractNQuartets(net::HybridNetwork, phylip::Phylip)\n\nExtracts quartet information from a topology. Also stores the observed site patter frequency information for each quartet.\n\nInput\n\nnet   HybridNetwork object obtained using the function readTopology\nphylip    Phylip object obtained using the function readPhylipFile!\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.extractQuartets","page":"PhyNEST.jl","title":"PhyNEST.extractQuartets","text":"extractQuartets(net::HybridNetwork)\n\nExtracts quartet information from a topology.\n\nInput\n\nnet   HybridNetwork object obtained using the function readTopology\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.printQuartets","page":"PhyNEST.jl","title":"PhyNEST.printQuartets","text":"printQuartets(net::HybridNetwork)\n\nVisualizes quartet information extracted from a topology in the form of DataFrame.\n\nInput\n\nnet HybridNetwork object obtained using the function readTopology\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.TrueSitePatternSymm","page":"PhyNEST.jl","title":"PhyNEST.TrueSitePatternSymm","text":"TrueSitePatternSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nTrueSitePatternSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the symmetric quartet tree, ((1,2),(3,4));.  Three speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default.  See manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nInput\n\nmyt1 Speciation time for the common ancestor of species 1 and 2\nmyt2 Speciation time for the common ancestor of species 3 and 4\nmyt3 Root age\ntheta Effective population size parameter\nalpha 4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.TrueSitePatternAsymm","page":"PhyNEST.jl","title":"PhyNEST.TrueSitePatternAsymm","text":"TrueSitePatternAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)\nTrueSitePatternAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\nComputes true site pattern probabilities for the asymmetric quartet tree, (1,(2,(3,4)));.  Three speciation times (or node ages) in coalescent unit and theta must be provided. Alhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default.  See manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nInput\n\nmyt1 Speciation time for the common ancestor of species 3 and 4\nmyt2 Speciation time for the common ancestor of species 2, 3 and 4\nmyt3 Root age\ntheta Effective population size parameter\nalpha 4/3 by default\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.PhyNE!","page":"PhyNEST.jl","title":"PhyNEST.PhyNE!","text":"PhyNE!(startT::HybridNetwork,inputFile::Phylip,outgroup::String)\n\nEstimates the network or tree to fit observed site pattern frequencies stored in a Phylip object, using composite likelihood. A level-1 network is assumed. The search begins from topolgoy startT, which can be a tree or a network (with less than hmax reticulation nodes). Usually, startT is expected to be a tree estiamted from data, but can also be randomly generated. The topology is rooted using the outgroup provided. This must be identical to the outgroup sequence provided in the  sequence alignment. By default, PhyNE! will search the network space using simulated annealing assuming  hmax=1.\n\nThere are many optional arguments (values in parenthesis are defaults):\n\nhillclimbing(=false) Select hill climbing search\nhmax(=1) Maximum number of reticulations in the final network\nnruns(=10) Number of independent runs\nnniStartT=false \nmaxcount=1000 Number of steps to terminate the search\nnfail=75 Number of consecutive failures to terminate the search\n\nSimulated annealing stuff:\n\nburninn=25\nk=10\nProbabilityOfNotSelectingNode=9.5e-45\ncons=0.9\nalph=0.8\n\nOptimization stuff:\n\nNumIter=1000\nparamprint=false\n\nMiscellaneous:\n\ntimestamp=false\nfilename=\"\"\ndisplay=falses\n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNEST.jl","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    Phylip\n    nquartets\n    Nquartets","category":"page"},{"location":"#PhyNEST.Phylip","page":"PhyNEST.jl","title":"PhyNEST.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type. Sequence alignment with the following attributes:\n\nfilename    The actual filename of the input sequence alignment\n\ntime    Time it took to intake the filename in seconds\n\nnumtaxa The number of taxa appears on the first line of phylip\n\nseqleng The sequence length appears on the first line of phylip\n\nnametaxa    Names of each individual in the order of appearance in phylip\n\ncounttaxa   Numbering each individual in the order of appearance in phylip\n\nallquartet  All combinations of quartets using nametaxa\n\nindex   Just convert allquartet elements into arbitrary index numbers\n\nspcounts    Arrays of 15 site pattern frequencies \n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.nquartets","page":"PhyNEST.jl","title":"PhyNEST.nquartets","text":"nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a quartet extracted from a topology with the following attributes:\n\nnumber  List of individuals quartets in the order of i,j,k,l\n\ntquartet    List of quartets in the order of i,j,k,l\n\nindexNET    Index number for each quartet in the order appears in quartnet. \n\nmspcountsNET    We can directly use this counts for likelihood calculation.\n\nmrca    List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\n\nntau    Unique number of the taus in the tau used in branchlengths\n\nbranchlength    Branch length from the mrca to tau=0 in the order of i and j (ij),ik,il,jk,jl,and kl\n\nmomestlength    Branch length identified for each quartet given the mspcounts using moment estimator - Will get filled later because theta is required\n\nmombl   Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\n\nsymtype Type of each quartet. It can be either symmetric or asymmetric. Asymmetric (type 0) can have four possible topologies and we define type 1 as (i,((j,k),l)); type 2 as ((i,(j,k)),l); type 3 as (i,(j,(k,l))); and type 4 as (((i,j),k),l).\n\nlogLik  It's just a likelihood for the quartet\n\n\n\n\n\n","category":"type"},{"location":"#PhyNEST.Nquartets","page":"PhyNEST.jl","title":"PhyNEST.Nquartets","text":"Nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a tree extracted from a topology with the following attributes:\n\nleafname    Name of the leaves as appears on the newick\n\nleafnumber  Assigned number for each leaf when reading in the newick\n\ngamma   Gamma as appears on the newick, although irrelevant for our computation\n\nnquartet    An array of all quartets in the given topology\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNEST.jl","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    extractNQuartets\n    printQuartets\n    TrueSitePatternSymm\n    TrueSitePatternAsymm\n    PhyNE!\n    Phylip\n    nquartets\n    Nquartets","category":"page"}]
}
