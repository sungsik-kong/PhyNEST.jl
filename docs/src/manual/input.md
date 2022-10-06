# Input

Currently, PhyNE takes a relaxed, sequential, PHYLIP-formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. For more information, see [here](https://en.wikipedia.org/wiki/PHYLIP). 

## Reading in DNA sequence data
The following tutorial will use an example DNA sequence alignment file, `n5h1_3k.txt`, available in `example` folder of the software. This alignment is 750000 base pairs (bp) long generated from 3000 simulated gene trees using a single hybridization scenario. Five sequences are in the alignment and named as `[1, 2, 3, 4, 5]`. 

After loading the package in Julia using the command `using PhyNE`, the folder where the alignment is stored must be specified. We set the location of the alignment `n5h1_3k.txt` as `path`. Then we can use the function `readPhylipFile!(path)` and specify the alignment. When the boolean option `showProgress` is set as true, PhyNE will visualize the progress of parsing the alignment. Here we set it as false for brevity.

```@example input
using PhyNE
path = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_3k.txt");
nothing # hide
```
```@repl input
data = readPhylipFile!(path, showProgress=false)
```

## Exporting the parsed site patterns frequencies in `.csv` file
```@repl input
data = readPhylipFile!(path, writecsv=true, showProgress=false)
```
```@repl input
csvpath = joinpath(dirname(pathof(PhyNE)), "..","example","sitePatternCounts_n5h1_3k.txt.csv");
df = readCSVFile(csvpath)
```

## Reading in `.ckp` file
Every time a sequence alignment is parsed, PhyNE creates a `.ckp` file that contains all information in the object `PHYPLIP`. PhyNE can parse a DNA alignment reasonably fast, however, it can take a while if the dataset is large. In case multiple network analysis using the same data is planned, one can bypass data parsing and use the automatically created `.ckp` file. 
```@repl input
ckppath = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_3k.txt.ckp");
ckpdata = readCheckPoint(ckppath)
```

