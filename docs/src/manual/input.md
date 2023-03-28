# Input

Current version of [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) takes a relaxed, sequential, [PHYLIP](https://en.wikipedia.org/wiki/PHYLIP) formatted DNA alignment file. 

A sample DNA alignment file `sample_n5h1.txt` is provided in the `/example` folder of the package [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) (or can be downloaded [here](https://github.com/sungsik-kong/PhyNE.jl/blob/main/example/n5h1_5k.txt)). The alignment contains five sequences `1, 2, 3, 4, and 5` and is 1,000,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): `(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

## Parsing DNA alignment data

Load [PhyNEST.jl](https://github.com/sungsik-kong/PhyNEST.jl) in Julia and move to the working directory that contains the alignment. The alignment is parsed with the function `readPhylip()`. The progress can be visualized with `showProgress` option that takes a boolean argument (set as `false` in this example for brevity). In the following command, we also set `checkpoint` option as `true` to create the `.ckp` file. Once the parsing the alignment is complete, a checkpoint file `n5h1_5k.txt.ckp` should appear in the directory. Note the `.ckp` file will have the same name as the alignment file. 

```@repl input
using PhyNEST
datapath = joinpath(dirname(pathof(PhyNEST)), "..","example","sample_n5h1.txt");
phydata = readPhylip(datapath, showProgress=false, checkpoint=true)
```

## Reading `.ckp` file
In case of multiple network analyses using the same alignment, a user can bypass the data parsing process by importing the `.ckp` file that is created after parsing the alignment for the first time. 

Reading a checkpoint file using the function `readCheckPoint()` as shown below.
```@repl input
ckppath = joinpath(dirname(pathof(PhyNEST)), "..","example","sample_n5h1.txt.ckp");
ckpdata = readCheckPoint(ckppath)
```