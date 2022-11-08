# Network estimation

## Getting ready for network estimation
To start the network analysis, (1) sequence alignment file is required and (2) starting topology and (3) an outgroup must be specified.

We use the `PHYLIP` alignment file in `/example` folder of PhyNEST package. The sequence alignment is parsed using function `readPhylipFile!` as shown [here](https://sungsik-kong.github.io/PhyNEST.jl/dev/manual/input/#Parsing-DNA-alignment-data). 

Then, we specify the starting tree (or network) stored in the `HybridNetwork` object using function [`readTopology()`](https://crsl4.github.io/PhyloNetworks.jl/latest/lib/public/#PhyloNetworks.readTopology). The starting tree can be either estimated from the data or randomly generated. We use a random tree topology.

```@julia netest
using PhyNEST
data = readPhylipFile!("n5h1_5k.txt")
startingtree = readTopology("(5,(4,(3,(2,1))));")
```

## Simulated Annealing
Function `PhyNE!()` executes the network analysis. By default, it uses simulated annealing algorithm to search the network space. `PhyNE!()` has a number of optional arguments (see [here](https://sungsik-kong.github.io/PhyNEST.jl/dev/#PhyNEST.PhyNE!)), but specifying the starting tree, data, and an outgroup and initiate the search as shown below: 
```@julia netest
netSA = PhyNE!(startingtree,data,"5")
```
To visualize the progress, the optional argument `display=true` can be used, which should print the following prompt, for example, at the initiation of the search.

    PhyNEST: Phylogenetic Network Estimation using SiTe patterns
    Analysis start: 2022-11-07 at 21:21:22
    Input filename: n5h1_5k.txt
    Number of sequences: 5 
    Sequence length: 2500000
    Starting Topology: (5,(4,(3,(2,1))));
    Outgroup specified for rooting: 5
    Number of maximum reticulation(s): 1
    The maximum number of iterations for each optimization: 1000
    Search algorithm selected: Simulated Annealing
    The maximum number of steps during search: 100000
    Alpha: 0.8; Cons: 0.9

    Initiating 5 iterations...

One of two output files `PhyNEST.sa.log` records all log throughout the search. For each iteration (or independent run) using simulated annealing search, PhyNEST records *k* best networks searched along with other relevant information as shown below:

    (1/5) Searching for the best network using the simulated annealing algorithm...
    Starting topology modified to (5,(4,(3,(2,1))));
    Running 25 runs of burn-in...Complete 
    (Cooling schedule)U=299136.023507182
    (Cooling schedule)beta=0.8999957157391261
    Rank   Composite Likelihood    Network
    1	2.88360994739e6         (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);
    2	2.88360994739e6         (5:10.143246111976964,(4:5.128207274589672,((3:1.7832993566222581,(2)#H6:::0.5947475713189532):1.3303497322531534,(#H6:::0.40525242868104683,1:1.1961994949262897):1.9174495939491218):2.01455818571426):5.0150388373872925);
    3	2.88362151446e6         (5:10.694491971616939,(4:5.43212278921928,((3:2.485548434584157,(#H6:::0.32538307301269104,2:6.112907457836508e-13):2.4855484345835457):1.2278094075673045,(1)#H6:::0.674616926987309):1.7187649470678181):5.262369182397659);
    4	2.88362151446e6         (5:10.694492301887275,(4:5.43212302522998,((3:2.485548536348796,(2:3.738582439844855e-12,(1)#H6:::0.32538305925513744):2.4855485363450573):1.2278094281091798,#H6:::0.6746169407448626):1.7187650607720037):5.262369276657295);
    5	2.88362693623e6         (5:10.725629067294516,(4:5.4494219219870095,(#H6:::0.6715226076537066,(1:2.508818645237626,(2:6.922689088873185e-13,(3)#H6:::0.32847739234629336):2.508818645236934):1.2147065868268707):1.7258966899225126):5.276207145307507);
    6	2.88362693623e6         (5:10.725629357953252,(4:5.449422058791485,((1:2.508818707391581,(#H6:::0.32847740524951685,2:2.954945769153218e-13):2.5088187073912858):1.2147066492672614,(3)#H6:::0.6715225947504831):1.7258967021326423):5.276207299161768);
    7	2.8837567506e6         (5:12.935996175051075,(#H6:::0.21325293589105276,((1:3.1629776123192412,(2:1.933312084245944,(3)#H6:::0.7867470641089472):1.2296655280732973):3.365071192479882,4:6.528048804799123):1.170591428480983):5.237355941770969);
    8	2.88378761042e6         (5:13.131954265122538,((3:3.2109280989424116,(2:1.8718749279421563,(1)#H6:::0.7233614942899573):1.3390531710002553):3.5626956885379393,(4:6.773623787480175,#H6:::0.27663850571004267):1.758593271006248e-13):6.358330477642187);
    9	2.88380850997e6         (5:13.030209792938253,((1:3.198422341899994,(2:1.8541586155663632,(3)#H6:::0.7273145740149164):1.3442637263336308):3.5191039964160598,(#H6:::0.2726854259850836,4:6.717526338316051):2.6645352591003757e-15):6.312683454622199);
    10	2.88380850997e6         (5:13.030209981281956,(#H6:::0.27268542338357454,(4:6.71752636287232,(1:3.198422435121352,(2:1.8541587632941887,(3)#H6:::0.7273145766164255):1.3442636718271634):3.519103927750968):0.0):6.312683618409636);
    Speciation times for some newicks may not have updated if estimates are weird (e.g., NaN).
    The search terminated at step 360 and at 50th consecutive failures and (Cooling schedule)ci=922.9788432411981.
    Summary of each move:
    Insertion of reticulation edge: 1 proposed, 1 accepted.
    Tail move of reticulation edge: 99 proposed, 37 accepted. 
    Head move of reticulation edge: 101 proposed, 11 accepted.
    Change the direction of reticulation edge: 82 proposed, 19 accepted.
    Deletion of reticulation edge: 0 proposed, 0 accepted.
    Nearest-neighbor interchange (NNI): 77 proposed, 9 accepted.
    On the current topology, 60 moves were made, including 10 unsuccessful moves.
    Terminated because it reached the maximum number of failures (current nfail=50).
    The best network found in this run: (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,(2)#H6:::0.4052526191450872):1.9174483823877784,(#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);
    -Log Composite Likelihood: 2.8836099473859877e6 

## Hill climbing
Network search using hill climbing algorithm can be initiated using an option `hillclimbing=true`.

```@julia netest
netHC = PhyNE!(startingtree,data,"5",hillclimbing=true)
```
Similar to the simulation annealing, it will create two output files named `PhyNe.hc.log` and `PhyNe.hc.out`





## Network visualization
At the end of the search, the output file with an extension `.out` is created in the working directory that should look something similar to below:

    The best network found from 5 runs using the simulated annealing algorithm
    MCL network: 
    (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,#H6:::0.4052526191450872):1.9174483823877784,((2)#H6:::0.5947473808549129,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);
        Dendroscope: 
    (5:10.143244413452512,(4:5.128206352203117,((1:1.1962000978367888,#H6):1.9174483823877784,((2)#H6,3:1.783298580712527):1.3303498995120402):2.0145578719785497):5.015038061249395);
        -Log Composite Likelihood: 2.8836099473859877e6.
    end

Two extended Newick strings represent the identical network topology, but formatted for visualization using the Julia package [`PhyloPlots`](https://github.com/cecileane/PhyloPlots.jl) (top) or [Dendroscope 3](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/dendroscope/)(bottom).