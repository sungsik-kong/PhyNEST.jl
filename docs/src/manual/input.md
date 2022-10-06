# Input

Currently, PhyNE takes a relaxed, sequential, PHYLIP-formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. For more information, see [here](https://en.wikipedia.org/wiki/PHYLIP). 

## Reading in DNA sequence data
The following tutorial will use an example DNA sequence alignment file, `n5h1_3k.txt`, available in `example` folder of the software. This alignment is 750000 base pairs (bp) long generated from 3000 simulated gene trees using a single hybridization scenario. Five sequences are in the alignment and named as `[1, 2, 3, 4, 5]`. 

After loading the package in Julia using the command `using PhyNE`, the folder where the alignment is stored must be specified. We set the location of the alignment `n5h1_3k.txt` as `path`. Then we can use the function `readPhylipFile!(path)` and specify the alignment. When the boolean option `showProgress` is set as true, PhyNE will visualize the progress of parsing the alignment. Here we set it as false for brevity.

```@example qcf
using PhyNE
path = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_3k.txt");
nothing # hide
```
```@repl qcf
data = readPhylipFile!(path, showProgress=false)
```

## Reading in `.ckp` file

## Exporting as `.csv` file
