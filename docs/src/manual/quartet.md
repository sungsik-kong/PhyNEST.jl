# Quartet

## Reading a network
An extended Newick string can be read by the function `readTopology` from `PhyloNetworks`. Here we will use a five tip network:

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

```@repl quartet
using PhyNE
network = readTopology("(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));")
```
## Extract quartet(s)

net=readTopology
```@repl quartet
printQuartets(extractQuartets(network))
```

## True site pattern probabilities for a quartet
### Symmaetric quartet
GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
net=readTopology
```@repl quartet
symqProb=GetTrueProbsSymm(1.0,2.0,5.0,0.01,4/3)
```
### Asymmetric quartet
function GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
```@repl quartet
asymqProb=GetTrueProbsAsymm(1.0,2.0,5.0,0.01,4/3)
```
## Simulate true site pattern frequencies
function simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)
```@repl quartet
simSPsym=simspcounts(0,1.0,2.0,5.0,0.01,4/3,1000000)
simSPasym=simspcounts(3,1.0,2.0,5.0,0.01,4/3,1000000)
```
## Method-of-moment estimator of branch lengths
function momentEstimat(type::Integer,spcounts::Array,theta::Float64)
```@repl quartet
momEstsym=momentEstimat(0,simSPsym,0.01)
momEstasym=momentEstimat(3,simSPasym,0.01)
```

## Estimating theta
function startTheta(q::Array{Nquartets, 1},net::HybridNetwork; lbound=0.00001::Float64,factor=2.0::Float64)
