var documenterSearchIndex = {"docs":
[{"location":"manual/input/#Input","page":"Input","title":"Input","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Current version of PhyNEST.jl takes a relaxed, sequential, PHYLIP formatted DNA alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"A sample DNA alignment file sample_n5h1.txt is provided in the /example folder of the package PhyNEST.jl (or can be downloaded here). The alignment contains five sequences 1, 2, 3, 4, and 5 and is 1,000,000 base pairs (bp) long. The true relationship of the five taxa is written in the extended Newick format as (branch lengths not shown): (5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));.","category":"page"},{"location":"manual/input/#Parsing-DNA-alignment-data","page":"Input","title":"Parsing DNA alignment data","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"Load PhyNEST.jl in Julia and move to the working directory that contains the alignment. The alignment is parsed with the function readPhylip(). The progress can be visualized with showProgress option that takes a boolean argument (set as false in this example for brevity). In the following command, we also set checkpoint option as true to create the .ckp file. Once the parsing the alignment is complete, a checkpoint file n5h1_5k.txt.ckp should appear in the directory. Note the .ckp file will have the same name as the alignment file. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"using PhyNEST\ndatapath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"sample_n5h1.txt\");\nphydata = readPhylip(datapath, showProgress=false, checkpoint=true)","category":"page"},{"location":"manual/input/#Reading-.ckp-file","page":"Input","title":"Reading .ckp file","text":"","category":"section"},{"location":"manual/input/","page":"Input","title":"Input","text":"In case of multiple network analyses using the same alignment, a user can bypass the data parsing process by importing the .ckp file that is created after parsing the alignment for the first time. ","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"Reading a checkpoint file using the function readCheckPoint() as shown below.","category":"page"},{"location":"manual/input/","page":"Input","title":"Input","text":"ckppath = joinpath(dirname(pathof(PhyNEST)), \"..\",\"example\",\"sample_n5h1.txt.ckp\");\nckpdata = readCheckPoint(ckppath)","category":"page"},{"location":"manual/installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"manual/installation/#Julia-installation","page":"Installation","title":"Julia installation","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Julia aims to create an unprecedented combination of ease-of-use, power, and efficiency in a single programming language. Current release of Julia is available here. Julia can be credited by citing the following paper:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Jeff Bezanson, Alan Edelman, Stefan Karpinski, Viral B. Shah. (2017) Julia: A Fresh Approach to Numerical Computing. SIAM Review, 59: 65–98. doi: 10.1137/141000671.","category":"page"},{"location":"manual/installation/#Installation-of-the-package-PhyNEST","page":"Installation","title":"Installation of the package PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To use PhyNEST.jl, Julia >= v1.7 is recommended. This package was developed in Julia 1.7.2, and has been tested for Julia >= v1.7 in OSX distributions.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To install PhyNEST.jl run the following command in Julia prompt:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using Pkg\nPkg.add(\"PhyNEST\")","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"or","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] add PhyNEST","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Either command installs the package. To check the current version of installed PhyNEST.jl, run the following command:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"] status PhyNEST","category":"page"},{"location":"manual/installation/#Update-PhyNEST","page":"Installation","title":"Update PhyNEST","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"To update PhyNEST.jl to the most recent version, run:","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Pkg.update(\"PhyNEST\")","category":"page"},{"location":"manual/installation/#Test-example","page":"Installation","title":"Test example","text":"","category":"section"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"Once the installation is complete, load PhyNEST.jl and run function greet() to see if the package has been added to the local machine. A greet message with current time should appear.","category":"page"},{"location":"manual/installation/","page":"Installation","title":"Installation","text":"using PhyNEST\ngreet()","category":"page"},{"location":"#PhyNEST.jl","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Documentation for a Julia package PhyNEST.jl: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"","category":"page"},{"location":"#References","page":"PhyNEST.jl","title":"References","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.","category":"page"},{"location":"#Getting-help","page":"PhyNEST.jl","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Manual","page":"PhyNEST.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"PhyNEST.jl","title":"PhyNEST.jl","text":"Pages = [\n    \"manual/installation.md\",\n    \"manual/input.md\",\n    ]","category":"page"}]
}
