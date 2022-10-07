# Quartet

## Reading a network
An extended Newick string can be read by the function `readTopology` from `PhyloNetworks`. Here we will use a five tip network:

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

```@repl quartet
network = readTopology("(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));")
```
## Extract quartet(s)
net=readTopology
q=extractNQuartets2(net)
print(q)

## True site pattern probabilities for a quartet
### Symmaetric quartet
GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
### Asymmetric quartet
function GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)

## Simulate true site pattern frequencies
function simspcounts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)

## Method-of-moment estimator of branch lengths
function momentEstimat(type::Integer,spcounts::Array,theta::Float64)


## Estimating theta
