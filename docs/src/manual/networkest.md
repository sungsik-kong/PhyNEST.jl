# Network estimation
To start the network analysis, we need sequence alignment file, starting topology (that can be either randomly generated or estimated from the data using PhyNEST or other method), and an outgroup. In this example, we are going to use the Phylip alignment in the `/example` folder of PhyNEST package. First, we parse the alignment using `readPhylipFile!` function followed by reading in a random topology using `readTopology` function (see below).

```@julia netest
using PhyNEST
data = readPhylipFile!("n5h1_5k.txt")
startingtree = readTopology("(5,(4,(3,(2,1))));")
```

## Simulated Annealing
Function `PhyNE!()` executes the network analysis. By default, it uses simulated annealing alogrithm to search the network space. While there are a number of options, the search can be initiated with the starting topology and data with an outgroup specifed. Below command can be used, for example. 

```@julia netest
netSA = PhyNE!(startingtree,data,"outgroup")
```
To visualize the progress, an option `display=true` can be used, and the following prompt should be printed on the screen that summarizes some information of the data and the searching condition.

    PhyNE: Estimating Maximum Pseudolikelihood Phylogenetic Network
    ╔═══╦╗─────╔═╗─╔╗
    ║╔═╗║║─────║║╚╗║║
    ║╚═╝║╚═╦╗─╔╣╔╗╚╝╠══╗
    ║╔══╣╔╗║║─║║║╚╗║║║═╣
    ║║──║║║║╚═╝║║─║║║║═╣
    ╚╝──╚╝╚╩═╗╔╩╝─╚═╩══╝
    ───────╔═╝║ 
    ───────╚══╝   
    Analysis start: 2022-11-02 at 17:18:20
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

Once the analysis is complete, two output files (named `PhyNe.sa.log` and `PhyNe.sa.out` by default) are stored. The `.log` file contains all output from each independent run. The `.out` file contains the network topology found with the best composite likelihood, written in formats readable in the package PhyloNetworks and Dendroscope, for visualization.


## Hill climbing
Network search using hill climbing algorithm can be initiated using an option `hillclimbing=true`.

```@julia netest
netHC = PhyNE!(startingtree,data,"5",hillclimbing=true)
```
Similar to the simulation annealing, it will create two output files named `PhyNe.hc.log` and `PhyNe.hc.out`