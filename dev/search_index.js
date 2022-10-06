var documenterSearchIndex = {"docs":
[{"location":"manual/input/#Input","page":"Input","title":"Input","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Currently, PhyNE takes a relaxed, sequential, PHYLIP-formatted DNA alignment file. Extensions can be either .phy or .txt. For more information, see here. ","category":"page"},{"location":"manual/input/#Reading-in-DNA-sequence-data","page":"Input","title":"Reading in DNA sequence data","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"The following tutorial will use an example DNA sequence alignment file, n5h1_3k.txt, available in example folder of the software. This alignment is 750000 base pairs (bp) long generated from 3000 simulated gene trees using a single hybridization scenario. Five sequences are in the alignment and named as [1, 2, 3, 4, 5]. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"After loading the package in Julia using the command using PhyNE, the folder where the alignment is stored must be specified. We set the location of the alignment n5h1_3k.txt as path. Then we can use the function readPhylipFile!(path) and specify the alignment. When the boolean option showProgress is set as true, PhyNE will visualize the progress of parsing the alignment. Here we set it as false for brevity.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"using PhyNE\npath = joinpath(dirname(pathof(PhyNE)), \"..\",\"example\",\"n5h1_3k.txt\");\nnothing # hide","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"data = readPhylipFile!(path, showProgress=false)","category":"page"},{"location":"manual/input/#Exporting-the-parsed-site-patterns-frequencies-in-.csv-file","page":"Input","title":"Exporting the parsed site patterns frequencies in .csv file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"data = readPhylipFile!(path, writecsv=true, showProgress=false)","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"csvpath = joinpath(dirname(pathof(PhyNE)), \"..\",\"example\",\"sitePatternCounts_n5h1_3k.txt.csv\");\ndf = readCSVFile(csvpath)","category":"page"},{"location":"manual/input/#Reading-in-.ckp-file","page":"Input","title":"Reading in .ckp file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Every time a sequence alignment is parsed, PhyNE creates a .ckp file that contains all information in the object PHYPLIP. PhyNE can parse a DNA alignment reasonably fast, however, it can take a while if the dataset is large. In case multiple network analysis using the same data is planned, one can bypass data parsing and use the automatically created .ckp file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"ckppath = joinpath(dirname(pathof(PhyNE)), \"..\",\"example\",\"n5h1_3k.txt.ckp\");\nckpdata = readCheckPoint(ckppath)","category":"page"},{"location":"manual/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Manual TBA","category":"page"},{"location":"manual/networkest/#Network-estimation","page":"Network estimation","title":"Network estimation","text":"","category":"section"},{"location":"manual/networkest/#Optimizing-parameters-for-a-given-network","page":"Network estimation","title":"Optimizing parameters for a given network","text":"","category":"section"},{"location":"manual/networkest/#Network-analysis","page":"Network estimation","title":"Network analysis","text":"","category":"section"},{"location":"manual/networkest/#Hill-climbing","page":"Network estimation","title":"Hill climbing","text":"","category":"section"},{"location":"manual/networkest/#Simulated-Annealing","page":"Network estimation","title":"Simulated Annealing","text":"","category":"section"},{"location":"manual/networkest/#Visualizing-network","page":"Network estimation","title":"Visualizing network","text":"","category":"section"},{"location":"manual/others/#Others","page":"Others","title":"Others","text":"","category":"section"},{"location":"manual/others/#HyDe","page":"Others","title":"HyDe","text":"","category":"section"},{"location":"manual/others/#The-D-statistics","page":"Others","title":"The D-statistics","text":"","category":"section"},{"location":"manual/others/","page":"Others","title":"Others","text":"Manual TBA","category":"page"},{"location":"manual/quartet/#Quartet","page":"Quartet","title":"Quartet","text":"","category":"section"},{"location":"manual/quartet/#Reading-a-network","page":"Quartet","title":"Reading a network","text":"","category":"section"},{"location":"manual/quartet/#Extract-quartet(s)","page":"Quartet","title":"Extract quartet(s)","text":"","category":"section"},{"location":"manual/quartet/#True-site-pattern-probabilities-for-a-quartet","page":"Quartet","title":"True site pattern probabilities for a quartet","text":"","category":"section"},{"location":"manual/quartet/#Symmaetric-quartet","page":"Quartet","title":"Symmaetric quartet","text":"","category":"section"},{"location":"manual/quartet/#Asymmetric-quartet","page":"Quartet","title":"Asymmetric quartet","text":"","category":"section"},{"location":"manual/quartet/#Simulate-true-site-pattern-frequencies","page":"Quartet","title":"Simulate true site pattern frequencies","text":"","category":"section"},{"location":"manual/quartet/#Method-of-moment-estimator-of-branch-lengths","page":"Quartet","title":"Method-of-moment estimator of branch lengths","text":"","category":"section"},{"location":"manual/quartet/#Estimating-theta","page":"Quartet","title":"Estimating theta","text":"","category":"section"},{"location":"#PhyNE.jl","page":"PhyNE.jl","title":"PhyNE.jl","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"Documentation for a julia package PhyNE.jl: Phylogenetic Network Estimation using composite likelihood directly from the sequence data using site pattern frequency. ","category":"page"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"","category":"page"},{"location":"#References","page":"PhyNE.jl","title":"References","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"TBA","category":"page"},{"location":"#Getting-help","page":"PhyNE.jl","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"Please use google group to report bugs or make suggestions.","category":"page"},{"location":"#Manual","page":"PhyNE.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"Pages = [\n    \"manual/installation.md\",\n    \"manual/input.md\",\n    \"manual/quartet.md\",\n    \"manual/networkest.md\",\n    \"manual/others.md\",\n    ]","category":"page"},{"location":"#Functions","page":"PhyNE.jl","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    mrca\n    nature\n    extractNQuartets\n    printQuarts\n    GetTrueProbsSymm\n    GetTrueProbsAsymm\n    GetTrueProbsNetTypes\n    simspcounts\n    momentEstimat\n    startTheta\n    Optimization\n    hillClimb\n    maxburnin\n    simAnneal\n\n    PhyNe!\n\n    Dstat\n    HyDe\n    LRT","category":"page"},{"location":"#PhyNE.readCheckPoint","page":"PhyNE.jl","title":"PhyNE.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file. gegeege\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.writeSitePatternCounts","page":"PhyNE.jl","title":"PhyNE.writeSitePatternCounts","text":"writeSitePatternCounts(p::Phylip,write::Bool,inputfile::AbstractString)\nwriteSitePatternCounts(p::Phylip,inputfile::AbstractString)\nwriteSitePatternCounts(inputfile::AbstractString)\nwriteSitePatternCounts(p::Phylip)\n\nWrites the observed site pattern frequencies from the phylip file as .csv file in the working  directory and name it as sitePatternCounts_inputfile.csv.  Although it requires an attribute inputfile, it does not actually read the file as long as the Phylip object in memory is provided. (i.e., any AbstractString can be provided here that  is used as a suffix of the ) However, if only the input file name is provided,  it runs readPhylipFile! first then write csv. file.\n\nInput\n\np A Type Phylip object\nwrite  Boolean variable to export site pattern frequences in .csv\ninputfile Name of phylip file as a AbstractString\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.nature","page":"PhyNE.jl","title":"PhyNE.nature","text":"nature(Node::Int64,net::HybridNetwork)\n\nFunction to determine the nature (root, internal tree node, leaf, hybrid node) of node given the node number and Hybridnetwork.\n\nInput\n\nNode       An integer node number\n\nnet        A tree/network in HybridNetwork\n\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.extractNQuartets","page":"PhyNE.jl","title":"PhyNE.extractNQuartets","text":"extractNQuartets(net::HybridNetwork,p::Phylip)\n\ngege\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.printQuarts","page":"PhyNE.jl","title":"PhyNE.printQuarts","text":"printQuarts\n\ngegejaga\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.GetTrueProbsSymm","page":"PhyNE.jl","title":"PhyNE.GetTrueProbsSymm","text":"GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\ngege\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.GetTrueProbsAsymm","page":"PhyNE.jl","title":"PhyNE.GetTrueProbsAsymm","text":"GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\ngege\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.GetTrueProbsNetTypes","page":"PhyNE.jl","title":"PhyNE.GetTrueProbsNetTypes","text":"GetTrueProbsNetTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)\n\ngege\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.simspcounts","page":"PhyNE.jl","title":"PhyNE.simspcounts","text":"simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)\n\nThis function simulates the 15 site pattern probabilities for a quartet.\n\n- `type`   specify the type of a quartet. \n- `myt1`   Branch length\n- `myt2`   Branch length\n- `myt3`   Branch length from the root\n- `theta`  Theta value (default=0.01)\n- `alpha`  Alpha value (default=4/3)\n- `length` Sequence lengths in integers (default=1000000)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.Optimization","page":"PhyNE.jl","title":"PhyNE.Optimization","text":"Optimization\n\nhehe\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.hillClimb","page":"PhyNE.jl","title":"PhyNE.hillClimb","text":"function hillClimb(startT::HybridNetwork,p::Phylip,hmax::Integer,outgroup::String,maxcount::Integer,\nnfail::Integer,NumIter::Integer,paramprint::Bool,logfile::IO,display::Bool)\n\nExecutes Hill-Climbing Algorithm.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNE.PhyNe!","page":"PhyNE.jl","title":"PhyNE.PhyNe!","text":"PhyNe!(startT::HybridNetwork,inputFile::Phylip,outgroup::String)\n\nEstimates the network or tree to fit observed site pattern frequencies stored in a Phylip object, using composite likelihood. A level-1 network is assumed. The search begins from topolgoy startT, which can be a tree or a network (with less than hmax reticulation nodes). Usually, startT is expected to be a tree estiamted from data, but can also be randomly generated. The topology is rooted using the outgroup provided. This must be identical to the outgroup sequence provided in the  sequence alignment. By default, PhyNe! will search the network space using simulated annealing assuming  hmax=1.\n\nThere are many optional arguments (values in parenthesis are defaults):\n\nhillclimbing(=false) Select hill climbing search\nhmax(=1) Maximum number of reticulations in the final network\nnruns(=10) Number of independent runs\nnniStartT=false \nmaxcount=1000 Number of steps to terminate the search\nnfail=75 Number of consecutive failures to terminate the search\n\nSimulated annealing stuff:\n\nburninn=25\nk=10\nProbabilityOfNotSelectingNode=9.5e-45\ncons=0.9\nalph=0.8\n\nOptimization stuff:\n\nNumIter=1000\nparamprint=false\n\nMiscellaneous:\n\ntimestamp=false\nfilename=\"\"\ndisplay=falses\n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNE.jl","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"    Phylip\n    nquartets\n    Nquartets","category":"page"},{"location":"#PhyNE.Phylip","page":"PhyNE.jl","title":"PhyNE.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type. Sequence alignment with the following attributes:\n\nfilename    The actual filename of the input sequence alignment\n\ntime    Time it took to intake the filename in seconds\n\nnumtaxa The number of taxa appears on the first line of phylip\n\nseqleng The sequence length appears on the first line of phylip\n\nnametaxa    Names of each individual in the order of appearance in phylip\n\ncounttaxa   Numbering each individual in the order of appearance in phylip\n\nallquartet  All combinations of quartets using nametaxa\n\nindex   Just convert allquartet elements into arbitrary index numbers\n\nspcounts    Arrays of 15 site pattern frequencies \n\n\n\n\n\n","category":"type"},{"location":"#PhyNE.nquartets","page":"PhyNE.jl","title":"PhyNE.nquartets","text":"nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a quartet extracted from a topology with the following attributes:\n\nnumber  List of individuals quartets in the order of i,j,k,l\n\ntquartet    List of quartets in the order of i,j,k,l\n\nindexNET    Index number for each quartet in the order appears in quartnet. \n\nmspcountsNET    We can directly use this counts for likelihood calculation.\n\nmrca    List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\n\nntau    Unique number of the taus in the tau used in branchlengths\n\nbranchlength    Branch length from the mrca to tau=0 in the order of i and j (ij),ik,il,jk,jl,and kl\n\nmomestlength    Branch length identified for each quartet given the mspcounts using moment estimator - Will get filled later because theta is required\n\nmombl   Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\n\nsymtype Type of each quartet. It can be either symmetric or asymmetric. Asymmetric (type 0) can have four possible topologies and we define type 1 as (i,((j,k),l)); type 2 as ((i,(j,k)),l); type 3 as (i,(j,(k,l))); and type 4 as (((i,j),k),l).\n\nlogLik  It's just a likelihood for the quartet\n\n\n\n\n\n","category":"type"},{"location":"#PhyNE.Nquartets","page":"PhyNE.jl","title":"PhyNE.Nquartets","text":"Nquartets\n\nSubtype of abstract NQuartet type. Stores relevant information of a tree extracted from a topology with the following attributes:\n\nleafname    Name of the leaves as appears on the newick\n\nleafnumber  Assigned number for each leaf when reading in the newick\n\ngamma   Gamma as appears on the newick, although irrelevant for our computation\n\nnquartet    An array of all quartets in the given topology\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNE.jl","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNE.jl","title":"PhyNE.jl","text":"    readPhylipFile!\n    readCheckPoint\n    writeSitePatternCounts\n    extractNQuartets\n    printQuarts\n    Phylip\n    mrca\n    nature\n    nquartets\n    Nquartets\n    PhyNe!","category":"page"}]
}
