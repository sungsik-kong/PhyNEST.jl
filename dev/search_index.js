var documenterSearchIndex = {"docs":
[{"location":"#PhyNEST","page":"PhyNEST","title":"PhyNEST","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Documentation for a Julia package PhyNEST: Phylogenetic Network Estimation using SiTe patterns.","category":"page"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"","category":"page"},{"location":"#Tutorials","page":"PhyNEST","title":"Tutorials","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"A step-by-step wiki tutorial is available, that has been done for the 2023 Botany workshop.","category":"page"},{"location":"#Reference","page":"PhyNEST","title":"Reference","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood. Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.","category":"page"},{"location":"#Getting-help","page":"PhyNEST","title":"Getting help","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"Please use google group to report bugs or post questions and/or suggestions.","category":"page"},{"location":"#Functions","page":"PhyNEST","title":"Functions","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    PhyNEST.greet\n    PhyNEST.readPhylip\n    PhyNEST.readPhylipFile\n    PhyNEST.PhylipFileInfo\n    PhyNEST.getUniqueQuartets\n    PhyNEST.sitePatternCounts\n    PhyNEST.spRearrange\n    PhyNEST.show_sp\n    PhyNEST.sitePatternsToDF\n    PhyNEST.write_sp\n    PhyNEST.storeCheckPoint\n    PhyNEST.readCheckPoint\n\n    PhyNEST.GetTrueProbsSymm\n    PhyNEST.GetTrueProbsAsymm\n\n    PhyNEST.Dstat\n    PhyNEST.showall","category":"page"},{"location":"#PhyNEST.greet","page":"PhyNEST","title":"PhyNEST.greet","text":"greet()\n\nDisplays a greeting with citation information. No argument needed. \n\nExample\n\njulia> greet()\nThank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.\nPlease report bugs or make suggestions to https://groups.google.com/g/phynest-users.\nIf you conduct an analysis using PhyNEST, please cite:\nSungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood.\nPreprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readPhylip","page":"PhyNEST","title":"PhyNEST.readPhylip","text":"readPhylip(inputfile::AbstractString; \n            writecsv=false::Bool,\n            csvname=\"\"::AbstractString,\n            showProgress=true::Bool,\n            tdigits=3::Integer,\n            checkpoint=false::Bool)\n\nImport, read, and parse the input phylip file. File name must be specified as string. By default, progress bar is displayed, and .csv and .ckp files are NOT produced.  See optional arguments below and modify if needed.\n\nMandatory argument\n\ninputfile     Name of the phylip file as a string\n\nOptional arguments\n\nwritecsv       (default=false) A boolean arguement to allow writing site pattern frequencies in .csv file\ncsvname        (default=sitePatternCounts_inputfile.csv) A string that will be name of the .csv file \nshowProgress   (default=true) A boolean argument for visualizing the progress of alignment parsing\ntdigits        (default=3) Number of decimal points to record timme taken to parse the file in seconds\ncheckpoint     (default=false) A boolean argument to store the Phylip object as a .ckp file. (Warning: .ckp file can be large)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readPhylipFile","page":"PhyNEST","title":"PhyNEST.readPhylipFile","text":"readPhylipFile(inputfile::AbstractString,\n                writecsv::Bool,\n                csvname::AbstractString,\n                showProgress::Bool)\n\nExceuted while running the function readPhylip.  Fills in the attributes in the Phylip object, except for the Phylip.time. \n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.PhylipFileInfo","page":"PhyNEST","title":"PhyNEST.PhylipFileInfo","text":"PhylipFileInfo(inputfile::AbstractString, \n                p::Phylip, \n                showProgress::Bool)\n\nFunction that shortens?summarizes? the sequence alignment into two matrices: \n\nThe one with unique sites and \nAnother with how many times each column in the previous matrix occurs throughout the alignment.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.getUniqueQuartets","page":"PhyNEST","title":"PhyNEST.getUniqueQuartets","text":"getUniqueQuartets(p::Phylip)\n\nGetting some combinations for four sequences (quartets). This will result in n! quartets where n is the number of sequences. I believe a better way to do this exists.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sitePatternCounts","page":"PhyNEST","title":"PhyNEST.sitePatternCounts","text":"sitePatternCounts(p::Phylip,\n                ppbase::Array,\n                counts::Array)\n\nComputes observed site pattern frequencies from UniqueBase and BaseCounts obtained from the function PhylipFileInfo.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.spRearrange","page":"PhyNEST","title":"PhyNEST.spRearrange","text":"spRearrange(p::Phylip)\n\nRearranges the elemetns of each unique quartet obtained using getUniqueQuartets, and also suffles the site pattern frequencies accordingly. This prevents 24 redundant computations, saving computation time.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.show_sp","page":"PhyNEST","title":"PhyNEST.show_sp","text":"show_sp(p::Phylip)\n\nPretty name for displaying all quartet and the corresponding site pattern frequencies on screen in the DataFrame format. It may result in excessively long table when there are many sequences in the input file.\n\nMandatory argument\n\np     Phylip object\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.sitePatternsToDF","page":"PhyNEST","title":"PhyNEST.sitePatternsToDF","text":"sitePatternsToDF(p::Phylip)\n\nExtracts the quartet and site pattern information from the Phyliip object, and reorganizes them in the DataFrame format.  This function does all hard work work for show_sp.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.write_sp","page":"PhyNEST","title":"PhyNEST.write_sp","text":"write_sp(p::Phylip; \n        csvname=\"PhyNEST_sp\"::AbstractString)\n\nPretty name for exporting all quartet and the corresponding site pattern frequencies in a .csv format.\n\nMandatory argument\n\np         Phylip object \n\nOptional argument\n\ncsvname   Filename for the .csv output can be given, otherwise we use PhyNEST_sp.csv by default.\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.storeCheckPoint","page":"PhyNEST","title":"PhyNEST.storeCheckPoint","text":"storeCheckPoint(p::Phylip)\n\nStores the phylip object into a .ckp file so a user does not have to repeat the phylip parsing again if working with the  same dataset next time. By default it creates a .ckp file that has the same filename as the input phylip file. File size can get pretty large... \n\nMandatory argument\n\np         Phylip object\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.readCheckPoint","page":"PhyNEST","title":"PhyNEST.readCheckPoint","text":"readCheckPoint(ckpfile::AbstractString)\n\nReads in .ckp file and creates a filled in phylip object.\n\nMandatory argument\n\nckpfile   Name of the checkpoint file\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsSymm","page":"PhyNEST","title":"PhyNEST.GetTrueProbsSymm","text":"GetTrueProbsSymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64)\nGetTrueProbsSymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64,\n                alpha::Float64)\n\nComputes true site pattern probabilities for the symmetric quartet tree: ((1,2),(3,4));. The fifteen quartet site pattern probabilities are returned  in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nThree speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified.  See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nMandatory arguments\n\nmyt1      Speciation time for the common ancestor of species 1 and 2 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\nOptional argument\n\nalpha (dafault=4/3)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.GetTrueProbsAsymm","page":"PhyNEST","title":"PhyNEST.GetTrueProbsAsymm","text":"GetTrueProbsAsymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64)\nGetTrueProbsAsymm(myt1::Float64,\n                myt2::Float64,\n                myt3::Float64,\n                theta::Float64,\n                alpha::Float64)\n\nComputes true site pattern probabilities for the asymmetric quartet tree: (1,(2,(3,4)));. The fifteen quartet site pattern probabilities are returned  in the order of:\n\nAAAA \nAAAB \nAABA \nAABB \nAABC \nABAA \nABAB \nABAC \nABBA \nBAAA \nABBC \nCABC \nBACA \nBCAA \nABCD\n\nThree speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified.  See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\n\nMandatory arguments\n\nmyt1      Speciation time for the common ancestor of species 3 and 4 in coalescent unit\nmyt2      Speciation time for the common ancestor of species 2 and (3,4) in coalescent unit\nmyt3      Root node age in coalescent unit\ntheta     Effective population size parameter\n\n##Optional argument\n\nalpha (dafault=4/3)\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.Dstat","page":"PhyNEST","title":"PhyNEST.Dstat","text":"Dstat(outgroup::String, p::Phylip;\n    pval=0.05::Float64, \n    displayall=false::Bool)\n\nConducts Patterson's D-statistic test. The result prints site pattern frequencies ABAB and ABBA used to compute  the D-statistic, Z-score, and the p-value for each quartet tested. Significance is marked with an asterisk. Function showall(df) can be subsequently used to show all rows.\n\nMandatory arguments\n\noutgroup     Name of the outgroup taxa\np   The Phylip object\n\nOptional arguments\n\npval       (default=0.05) Alpha level for significance\ndisplay_all (default=false) If set as true, the function print test results for every quartet. By default, it only prints those quartets where signficance was found.\n\nExample\n\njulia> p=readPhylip(\"sample_n4h1.txt\")\njulia> df=Dstat(\"4\",p)\nTip: if neccessary, use showall(df) function to see all the rows.\n2×10 DataFrame\n Row │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat     Zscore   pvalue   significance\n     │ String    String  String  String  Int64  Int64  Float64   Float64  Float64  String\n─────┼──────────────────────────────────────────────────────────────────────────────────────────\n   1 │ 4         3       2       1        1427   7852  0.692424  66.6995      0.0  *\n   2 │ 4         1       2       3        1427   7836  0.691892  66.5908      0.0  *\n\njulia> df=Dstat(\"4\",p,display_all=true)\nTip: if neccessary, use showall(df) function to see all the rows.\n6×10 DataFrame\nRow │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat        Zscore      pvalue    significance\n    │ String    String  String  String  Int64  Int64  Float64      Float64     Float64   String\n────┼─────────────────────────────────────────────────────────────────────────────────────────────────\n  1 │ 4         3       1       2        7852   1427  -0.692424    -66.6995    1.0\n  2 │ 4         3       2       1        1427   7852   0.692424     66.6995    0.0       *\n  3 │ 4         1       3       2        7836   1427  -0.691892    -66.5908    1.0\n  4 │ 4         1       2       3        1427   7836   0.691892     66.5908    0.0       *\n  5 │ 4         2       3       1        7836   7852   0.00101989    0.127743  0.449176\n  6 │ 4         2       1       3        7852   7836  -0.00101989   -0.127743  0.550824\n\n\n\n\n\n","category":"function"},{"location":"#PhyNEST.showall","page":"PhyNEST","title":"PhyNEST.showall","text":"showall(df::DataFrame)\n\nPrint all rows of the DataFrame object.    \n\n\n\n\n\n","category":"function"},{"location":"#Types","page":"PhyNEST","title":"Types","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    Phylip","category":"page"},{"location":"#PhyNEST.Phylip","page":"PhyNEST","title":"PhyNEST.Phylip","text":"Phylip\n\nSubtype of abstract INPUT type with the following attributes:\n\nfilename      Name of the input phylip alignment file\n\ntime          Time taken for parsing the input in seconds\n\nnumtaxa       Number of taxa given at the first line of the input file\n\nseqleng       Sequence length given at the first line of the input file\n\nnametaxa      Sequence names given in the input file\n\ncounttaxa     Unique integer identifier given to each individual in the order of appearance in the input file\n\nallquartet    All combinations of quartets for counttaxa\n\nindex         Just convert allquartet elements into arbitrary index numbers\n\nspcounts      Arrays of 15 site pattern frequencies for each quartet in allquartet\n\n\n\n\n\n","category":"type"},{"location":"#Index","page":"PhyNEST","title":"Index","text":"","category":"section"},{"location":"","page":"PhyNEST","title":"PhyNEST","text":"    ","category":"page"}]
}
