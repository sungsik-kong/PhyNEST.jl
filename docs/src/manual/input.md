# Input

Current version of PhyNE takes a relaxed, sequential, PHYLIP formatted DNA alignment file. Extensions can be either `.phy` or `.txt`. For more information on PHYLIP, see [here](https://en.wikipedia.org/wiki/PHYLIP). After the PHYLIP file has been parsed, PhyNE creates a checkpoint file with the extension `.ckp`, which can be called later to prevent repetition of parsing the same data and save time. Moreover, the observed site pattern frequencies of every quartet from the alignment can be exported into a `.csv` file, for a user to explore for other purposes (e.g., manually computation of Patterson's D-statistic etc.). The following sections explain above procedures in detail.

## Parsing DNA alignment data
A sample DNA alignment file called `n5h1_5k.txt` is given in the `example` folder of the package PhyNE (or can be downloaded [here](https://github.com/sungsik-kong/PhyNE.jl/blob/main/example/n5h1_5k.txt)). This alignment is 750,000 base pairs (bp) long and contains five sequences, named as  `[1, 2, 3, 4, 5]`. The alignment is generated using *ms* and *seq-gen*. The true relationship of the five taxa is written in the extended Newick format as (branch lengths omittted): 

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

We begin with loading the package PhyNE in Julia.
```@example input
using PhyNE
```

Assuming that you are working in the directory that contains the alignment file, the alignment file can be read using the function `readPhylipFile!(inputfile)`. Since the variable `inputfile` is expected to be a string, the file name should be in the quotation marks. The path of the file can be set as shown below or simply the name of the alignment file can be specified (e.g., `readPhylipFile!("n5h1_5k.txt")`). A checkpoint file with an extension `.ckp` is automatically created in the working directory upon completion of parsing the alignment. When the boolean option `showProgress` is set as `true`, PhyNE will show the progress bar during the data parsing process. Here, we set it as false for brevity.

```@repl input
datapath = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_5k.txt");
phydata = readPhylipFile!(datapath, showProgress=false)
```

## Reading `.ckp` file
Every time a sequence alignment is parsed, PhyNE creates a checkpoint file with the extension `.ckp`. Note the `.ckp` file will have the same name as the alignment file. In case multiple network analyses using the same data is planned, a user can bypass the (potentially time-consuming) data parsing process by using the function `readCheckPoint(ckpfile)` instead of the function `readPhylipFile!(inputfile)`. 
```@repl input
ckppath = joinpath(dirname(pathof(PhyNE)), "..","example","n5h1_5k.txt.ckp");
ckpdata = readCheckPoint(ckppath)
```

## Exporting as `.csv` file
The observed quartet site pattern frequencies from the data can be exported as a `.csv` formatted file. This can be done when parsing the alignment using the function `readPhylipFile!` by setting the option `writecsv=true`. It will create a `.csv` file in the working directory with the filename `sitePatternCounts_$inputfile.csv`. If a user would like to save the `.csv` file in a different name, the file name can be specified in the option `csvname`. In the below example, we set the filename as "tutorial_n5h1" for the `.csv` file. The `.csv` file can be visualized using function `readCSVFile(csvfile)`.
```@repl input
phydata = readPhylipFile!(datapath, showProgress=false,
                        writecsv=true, csvname="tutorial_n5h1")
csvpath = joinpath(dirname(pathof(PhyNE)), "..","example","tutorial_n5h1.csv");
csvdf = readCSVFile(csvpath)
```