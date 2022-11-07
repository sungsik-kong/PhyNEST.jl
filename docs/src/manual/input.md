# Input

Current version of [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) takes a relaxed, sequential, [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) formatted DNA alignment file. [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) automatically creates `FILENAME.ckp` upon the completion of reading in the alignment, and has an option for a user to export the observed site pattern frequencies as `.csv` file. The following sections explain these in detail.

## Parsing DNA alignment data
A sample DNA alignment file `n5h1_5k.txt` is provided in the `/example` folder of the package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) (or can be downloaded [here](https://github.com/sungsik-kong/PhyNE.jl/blob/main/example/n5h1_5k.txt)). The alignment contains five sequences `1, 2, 3, 4, and 5` and is 2,500,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): 

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

Load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) in Julia and move to the working directory that contains the alignment. The alignment is loaded with the function `readPhylipFile!()`. The progress can be visualized with `showProgress` option which is set as `false` in this example for brevity. Once the parsing the alignment is complete, a checkpoint file `n5h1_5k.txt.ckp` should appear in the directory.
```@repl input
using PhyNEST
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
