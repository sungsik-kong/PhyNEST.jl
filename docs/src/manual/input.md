# Input

PhyNE takes a PHYLIP-formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. 

## Reading in DNA sequence data
A sample sequence files, `n5h1_3k.txt` and `n5h1_5k.txt`, are available in the `\example` folder of the software.
```julia
less("n5h1_3k.txt")
```

```@example qcf
using PhyNE
align = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_3k.txt");
nothing # hide
```

```@repl qcf
sequences = readPhylipFile!(align, showProgress=false)
```

## Reading in `.ckp` file

## Exporting it as `csv` file
