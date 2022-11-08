# Getting ready for network estimation
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













Once the analysis is complete, two output files (named `PhyNe.sa.log` and `PhyNe.sa.out` by default) are stored. The `.log` file contains all output from each independent run. The `.out` file contains the network topology found with the best composite likelihood, written in formats readable in the package PhyloNetworks and Dendroscope, for visualization.


## Hill climbing
Network search using hill climbing algorithm can be initiated using an option `hillclimbing=true`.

```@julia netest
netHC = PhyNE!(startingtree,data,"5",hillclimbing=true)
```
Similar to the simulation annealing, it will create two output files named `PhyNe.hc.log` and `PhyNe.hc.out`