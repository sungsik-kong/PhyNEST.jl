# Input

Current version of [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) takes a relaxed, sequential, [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) formatted DNA alignment file. [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) automatically creates `FILENAME.ckp` upon the completion of reading in the alignment, and has an option for a user to export the observed site pattern frequencies as `.csv` file. The following sections explain these in detail.

## Parsing DNA alignment data
A sample DNA alignment named `n5h1_5k.txt` is provided in the `/example` folder of the package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) (or can be downloaded [here](https://github.com/sungsik-kong/PhyNE.jl/blob/main/example/n5h1_5k.txt)). The alignment contains five sequences named `[1, 2, 3, 4, 5]` and is 2,500,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): 

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

First, load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) in Julia.
```@example input
using PhyNEST
```








Assuming that you are working in the directory that contains the alignment file, the alignment file can be read using the function `readPhylipFile!(inputfile)`. Since the variable `inputfile` is expected to be a string, the file name should be in the quotation marks. The path of the file can be set as shown below or simply the name of the alignment file can be specified (e.g., `readPhylipFile!("n5h1_5k.txt")`). A checkpoint file with an extension `.ckp` is automatically created in the working directory upon completion of parsing the alignment. When the boolean option `showProgress` is set as `true`, PhyNEST will show the progress bar during the data parsing process. Here, we set it as false for brevity.

```@repl input
datapath = joinpath(dirname(pathof(PhyNEST)), "..","example","n5h1_5k.txt");
phydata = readPhylipFile!(datapath, showProgress=false)
```

## Reading `.ckp` file
Every time a sequence alignment is parsed, PhyNEST creates a checkpoint file with the extension `.ckp`. Note the `.ckp` file will have the same name as the alignment file. In case multiple network analyses using the same data is planned, a user can bypass the (potentially time-consuming) data parsing process by using the function `readCheckPoint(ckpfile)` instead of the function `readPhylipFile!(inputfile)`. 
```@repl input
ckppath = joinpath(dirname(pathof(PhyNEST)), "..","example","n5h1_5k.txt.ckp");
ckpdata = readCheckPoint(ckppath)
```

## Exporting as `.csv` file
The observed quartet site pattern frequencies from the data can be exported as a `.csv` formatted file. This can be done when parsing the alignment using the function `readPhylipFile!` by setting the option `writecsv=true`. It will create a `.csv` file in the working directory with the filename `sitePatternCounts_$inputfile.csv`. If a user would like to save the `.csv` file in a different name, the file name can be specified in the option `csvname`. In the below example, we set the filename as "tutorial_n5h1" for the `.csv` file. The `.csv` file can be visualized using function `readCSVFile(csvfile)`. 

If a user forgot to set `writecsv=true`, the function `writeSitePatternCounts` can be used to export. See [here](https://sungsik-kong.github.io/PhyNE.jl/dev/#PhyNE.writeSitePatternCounts) for more information.
```@repl input
phydata = readPhylipFile!(datapath, showProgress=false,
                        writecsv=true, csvname="tutorial_n5h1")
csvpath = joinpath(dirname(pathof(PhyNEST)), "..","example","tutorial_n5h1.csv");
csvdf = readCSVFile(csvpath)
```
