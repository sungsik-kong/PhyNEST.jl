# Input

Current version of PhyNE takes a relaxed, sequential, PHYLIP formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. For more information on PHYLIP, see [here](https://en.wikipedia.org/wiki/PHYLIP). After the PHYLIP file has been parsed, PhyNE creates a checkpoint file with the extension `.ckp`, which can be called later to prevent repetition of parsing the same data and save time. Moreover, the observed site pattern frequencies of every quartet from the alignment can be exported into a `.csv` file, for a user to explore for other purposes (e.g., manually computation of Patterson's D-statistic etc.). The following sections explain above procedures in detail.

## Parsing DNA alignment data
A sample DNA alignment file called `n5h1_5k.txt` is given in the `example` folder of the package PhyNE. This alignment is 750,000 base pairs (bp) long and contains five sequences, named as  `[1, 2, 3, 4, 5]`. The alignment is generated using *ms* and *seq-gen*. The true relationship of the five taxa is written in the extended Newick format as: (5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));. Species 5 is the outgroup, and species 2 is the hybrid daughter of species 3 and 5.

First, we load the package PhyNE in Julia.
```@example input
using PhyNE
```

Then, 

The folder where the alignment is stored must be specified. We set the location of the alignment `n5h1_3k.txt` as `path`. Then we can use the function `readPhylipFile!(path)` and specify the alignment. When the boolean option `showProgress` is set as true, PhyNE will visualize the progress of parsing the alignment. Here we set it as false for brevity.

```@repl input
datapath = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_5k.txt");
phydata = readPhylipFile!(datapath, 
                        writecsv=true, csvname="tutorial_n5h1", 
                        showProgress=false)
```

## Reading `.ckp` file
Every time a sequence alignment is parsed, PhyNE creates a `.ckp` file that contains all information in the object `PHYPLIP`. PhyNE can parse a DNA alignment reasonably fast, however, it can take a while if the dataset is large. In case multiple network analysis using the same data is planned, one can bypass data parsing and use the automatically created `.ckp` file. 
```@repl input
ckppath = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_5k.txt.ckp");
ckpdata = readCheckPoint(ckppath)
```

## Exporting as `.csv` file
```@repl input
csvpath = joinpath(dirname(pathof(PhyNE)), "..","example","tutorial_n5h1.csv");
csvdf = readCSVFile(csvpath)
```