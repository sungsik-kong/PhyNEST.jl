# Input

PhyNE takes a PHYLIP-formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. 

## Tutorial data

A sample sequence files, `n5h1_3k.txt` and `n5h1_5k.txt`, are available in the `\example` folder of the software.
```julia
less("n5h1_3k.txt")
```

```@example qcf
using PhyNE
align = joinpath(dirname(pathof(PhyloNetworks)), "..","examples","n5h1_1.txt");
nothing # hide
```

```@repl qcf
sequences = readPhylipFile(align)
```