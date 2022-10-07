# Quartet

## Reading a network
An extended Newick string can be read by the function `readTopology` from `PhyloNetworks`. Here we will use a five tip network:

`(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));`.

```@repl input
network = readTopology("(5,(4,((3,(2)#H6:::0.6),(1,#H6:::0.4))));")
```
## Extract quartet(s)

