# Network estimation
After parsing the input sequence alignment into quartet site pattern frequency data, we can estimate the network using a starting tree/network, that can be either randomly generated or estimated from the data. 
```@julia netest
startingtree = readTopology("(5,(4,(3,(2,1))));")
data=readPhylipFile!("phylipfile.phy")
```

## Hill climbing
Network search using hill climbing algorithm can be initiated using:
```@julia netest
netHC = PhyNE!(startingtree,data,"outgroup",hillclimbing=true)
```

To visualize the progress, an option `display=true` can be used, and the following prompt should be printed on the screen that summarizes some information of the data and the searching condition.

    PhyNe: Estimating Maximum Pseudolikelihood Phylogenetic Network
    ╔═══╦╗─────╔═╗─╔╗
    ║╔═╗║║─────║║╚╗║║
    ║╚═╝║╚═╦╗─╔╣╔╗╚╝╠══╗
    ║╔══╣╔╗║║─║║║╚╗║║║═╣
    ║║──║║║║╚═╝║║─║║║║═╣
    ╚╝──╚╝╚╩═╗╔╩╝─╚═╩══╝
    ───────╔═╝║ 
    ───────╚══╝   
    Analysis start: 2022-10-07 at 16:08:02
    Input filename: n5h1_5k.txt
    Number of sequences: 5 
    Sequence length: 2500000
    Starting Topology: (5,(4,(3,(2,1))));
    Outgroup specified for rooting: 1
    Number of maximum reticulation(s): 1
    The maximum number of iterations for each optimization: 1000
    Search algorithm selected: Hill-climbing
    The maximum number of steps during search: 50


## Simulated Annealing
```@julia netest
netSA = PhyNE!(startingtree,data,"outgroup")
```

