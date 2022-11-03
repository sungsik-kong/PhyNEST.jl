# Quartet

## Reading a network
An extended Newick string can be read by the function `readTopology` from `PhyloNetworks`. See [Cardona et al.,(2008)](https://doi.org/10.1186/1471-2105-9-532) to see the first description of the extended Newick formatted network. Here we use a five tip network:

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

PhyNE uses the same the extended Newick format as in `PhyloNetworks`. Branch lengths can be specified using colon (`:`) as in a regular Newick string. For reticulation nodes, relevant information are specified in the order of ':length:bootstrap:gamma'. PhyNE does *not* require branch lengths to be specified, however, if gamma is specified, it will be set as 0.5 .

```@repl quartet
using PhyNEST
network = readTopology("(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));")
```
## Extract quartet(s)
PhyNE can extract quartets extracted from a tree or a network using the function `extractQuartets`.  
```@repl quartet
printQuartets(network)
```

## True site pattern probabilities for a quartet
See [Chifman and Kubatko (2015)](https://www.sciencedirect.com/science/article/pii/S0022519315001095?via%3Dihub) for more information.
### Symmetric quartet
True site pattern probabilities for the fifteen site pattern for a symmetric quartet can be computed using the function `TrueSitePatternSymm`. We need to specify at least four parameters, three values of node ages, $\tau$, in the order of [MRCA of species 1 and 2, MRCA of species 3 and 4, root age] and population size parameter $\theta$. In the following block, we computed the probabilities for $\tau_1=1.0$, $\tau_2=2.0$, $\tau_3=5.0$, and $\theta=0.0025$. The 15 patterns appear in the order of [AAAA,AAAB,AABA,AABB,AABC,ABAA,ABAB,ABAC,ABBA,BAAA,ABBC,CABC,BACA,BCAA,ABCD]. 

We assure that the computation is on the right track by multiplying the weight for each site pattern, and sum all of them together to get 1.0.

```@repl quartet
symqProb=TrueSitePatternSymm(1.0,2.0,5.0,0.0025)
weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
sum(symqProb.*weights)
```
### Asymmetric quartet
function GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
```@repl quartet
asymqProb=TrueSitePatternAsymm(1.0,2.0,5.0,0.0025,4/3) 
```
## Simulate true site pattern frequencies
function simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)
```@repl quartet
simSPsym=simspcounts(0,1.0,2.0,5.0,0.0025,4/3,1000000)
```
```@repl quartet
simSPasym=simspcounts(3,1.0,2.0,5.0,0.0025,4/3,1000000)
```
## Method-of-moment estimator of branch lengths
function momentEstimat(type::Integer,spcounts::Array,theta::Float64)
```@repl quartet
momEstsym=momentEstimat(0,simSPsym,0.0025)
```
```@repl quartet
momEstasym=momentEstimat(3,simSPasym,0.0025)
```

## Estimating theta
function startTheta(q::Array{Nquartets, 1},net::HybridNetwork; lbound=0.00001::Float64,factor=2.0::Float64)
```@repl quartet
datapath = joinpath(dirname(pathof(PhyNEST)), "..","example","n5h1_5k.txt");
phydata = readPhylipFile!(datapath, showProgress=false)
theta=startTheta(network,phydata)
```


