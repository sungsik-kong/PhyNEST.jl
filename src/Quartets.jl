#written by Sungsik Kong 2021-2022

# 2.Quartet

#Use PhyloNetworks.readTopology to take newick or extnd newick formatted topology in.
#x=readTopology("(5,(4,((1,(2)#H6:::0.5),(3,#H6:::0.5))));")

#See the parsed topology
#PhyloNetworks.printEverything(x)
function tauNum(x::Int64) return abs(x+1) end

function extractNQuartets1(net::HybridNetwork)
    leafname,leafnumber=getLeafInfo(net)
    q=Nquartets(leafname,leafnumber)#creates Nquartets type.
    nq=listAllQs(q)#creates nquartets type in the function. 

    for x in nq.nquartet
        mrcapush(net,x)#push in the most recent common ancestors for each pair in a quartet
        symType(x)#identifies the symmetricity of the quartet tree
        #pushtau(x,net)#get the taus for each tree and rearranges it in the order of lengths
        uniqtaus(x) #gives unique number for each tau by altering the tree node number
    end
    return nq
end

"""
    extractNQuartets(net::HybridNetwork,p::Phylip)

gege
"""
function extractNQuartets(net::HybridNetwork,p::Phylip)
    nq=Nquartets[]
    if net.numHybrids==0
        nqq=extractNQuartets1(net)
        moveSPcounts(nqq,p)
        push!(nq,nqq)
    else
        t=displayedTrees(net,0.0) 
        numDispTrees=length(t)
        numDispTrees>1 || error("We expect 2 or more displayed trees, but we got $(length(t)).")
        for n in 1:numDispTrees
            nqq=extractNQuartets1(t[n])
            moveSPcounts(nqq,p)
            push!(nq,nqq)
        end
    end
    return nq
end

function extractQuartets(net::HybridNetwork)
    nq=Nquartets[]
    if net.numHybrids==0
        nqq=extractNQuartets1(net)
        push!(nq,nqq)
    else
        t=displayedTrees(net,0.0) 
        numDispTrees=length(t)
        numDispTrees>1 || error("We expect 2 or more displayed trees, but we got $(length(t)).")
        gamma=gamArray(net)
        for n in 1:numDispTrees
            nqq=extractNQuartets1(t[n])
            nqq.gamma=gamma[n]
            push!(nq,nqq)
        end
    end
    return nq
end


## Identifying Quartet Types

#### Getting the set of most recent common ancestors (MRCAs) of a quartet

"""
    getLeafInfo

**Function getLeafInfo** extracts leafname and leafnumber stored HybridNetwork object. Then checks 
if all taxa names are unique to make sure there is no duplicate tips in the input network.

### Input
**`net`**       A tree/network in Type object PhyloNetworks.HybridNetwork\\

### Output examples

"""
function getLeafInfo(net::HybridNetwork)
    leafname=String[]
    leafnumber=Int64[]

    for n in net.leaf
        if (n.leaf) 
            push!(leafname,n.name)
            push!(leafnumber,n.number)
        else error("Node $(n.name) is not a leaf.")             
        end
    end
    
    #check for the presence of duplicate leaves
    length(unique(leafname))==length(leafname) || error("Some taxa appears more than once in the input network.")
    
    return leafname,leafnumber
end
#leafname,leafnumber=getLeafInfo(x)#try it

"""
    listAllQs

**Function listAllQs** gets combinations of four based on the number of leaves in HybridNetwork stored in Type object 
*Nquartets.leafnumber*. It is function *getUniqueQuartets* above, but takes in the order stored in the HybridNetworks. 
For each set of four individuals, a new type object *nquartets* is created and stored in *Nquartets.nquartet*. 
There should be (n choose 4) quartets and this is checked. Part of the function *extractNQuartets*.

### Input
- **`nq`**       A tree/network in Type object *Nquartets* with leafname,leafnumber filled in\\

### Output examples

"""
function listAllQs(Nq::Nquartets)
    count=0
    numind=length(Nq.leafnumber) 
    numind>0 || error("There is no leaf specified in the input network")
    
    for i in (1:numind)
        for j in (i+1:numind)
            for k in (j+1:numind)
                for l in (k+1:numind)
                    count+=1; q=nquartets(count);
                    push!(q.quartet, [Nq.leafnumber[i],Nq.leafnumber[j],Nq.leafnumber[k],Nq.leafnumber[l]]); 
                    push!(Nq.nquartet,q);
                end
            end
        end
    end

    #check if the number is binomial (S 4)
    length(Nq.nquartet)==binomial(numind, 4) || error("Number of quartets extracted seems weird.")

    return Nq
end

function GetChild(edge::Edge) edge.node[edge.isChild1 ? 1 : 2] end
function GetParent(edge::Edge) edge.node[edge.isChild1 ? 2 : 1] end

"""
    childnode

Considering that an edge *(u,v)*, wher *u* is the parent of *v* if there is a directed route from *u* to *v* [see manuscript], 
**Funtion childnode** returns the nearest *v* given *u* in HybridNetwork. This function utilizes two functions from PhyloNetworks 
*getParent* and *getChild*. While all non-root internal tree nodes should have two child nodes, this function returns whichever
appears first in the network. Similarly, **Funtion parentnode** returns *u* given *v* and HybridNetwork.

### Input
- **`parentnode`**       Designated number of the parent node in interest\\
- **`network`**          A tree/network in Type object *Nquartets* with leafname,leafnumber filled in\\

### Output examples

"""
function childnode(parentnode::Int64,network::HybridNetwork)
    for edge in network.edge
        if parentnode==GetParent(edge).number
           childnode=GetChild(edge).number
        return childnode
        end
    end
end

"""
    childnode

Considering that an edge *(u,v)*, wher *u* is the parent of *v* if there is a directed route from *u* to *v* [see manuscript], 
**Funtion childnode** returns the nearest *v* given *u* in HybridNetwork. This function utilizes two functions from PhyloNetworks 
*getParent* and *getChild*. While all non-root internal tree nodes should have two child nodes, this function returns whichever
appears first in the network. Similarly, **Funtion parentnode** returns *u* given *v* and HybridNetwork.

### Input
- **`childnode`**       Designated number of the child node in interest\\
- **`network`**         A tree/network in Type object *Nquartets* with leafname,leafnumber filled in\\

### Output examples

"""
function parentnode(childnode::Int64,network::HybridNetwork)
    for edge in network.edge
        if childnode==GetChild(edge).number
            parentnode=GetParent(edge).number
        return parentnode
        end
    end
end



function mrca(Node1::Int64,Node2::Int64,net::HybridNetwork)
    root=net.node[net.root].number
    root!==nothing || error("Cannot identify root in the input tree")
    route1=Int64[]; push!(route1,Node1)
    route2=Int64[]; push!(route2,Node2)

    while Node1!==root 
        Node1=parentnode(Node1,net)
        Node1!==nothing || error("No parent exists for $Node1 in the input tree.")
        typeof(Node1)==Int64 || error("Parent node of the node $Node1 is not a node.")
        push!(route1,Node1)
    end

    while Node2!==root
        Node2=parentnode(Node2,net)
        Node2!==nothing || error("No parent exists for $Node2 in the input tree.")
        typeof(Node2)==Int64 || error("Parent node of the node $Node2 is not a node.")
        push!(route2,Node2)
    end

    for step1 in 1:length(route1)
        for step2 in 1:length(route2)
            if route1[step1]==route2[step2]
                return route1[step1]
            end
        end
    end

end

function mrcapush(net::HybridNetwork,q::nquartets)

    for indiv in q.quartet
        length(indiv)==4 || error("Each quartet must contain four taxa, but there is a quartet that contains $(length(indiv)) taxa.")
        
        i=indiv[1]; j=indiv[2]; k=indiv[3]; l=indiv[4]
        for leaf in [i,j,k,l]
            for n in net.node
                if n.number==leaf n.leaf || error("mrcapush: Node i: $(leaf) is not a leaf.") end
            end
        end

        ij=mrca(i,j,net)
        ik=mrca(i,k,net)
        il=mrca(i,l,net)
        jk=mrca(j,k,net)
        jl=mrca(j,l,net)
        kl=mrca(k,l,net)

        push!(q.mrca,[ij,ik,il,jk,jl,kl])
    end

end

function symType(q::nquartets)
    polytomy1=[0,1,2,3,4] #=(1,2,3,4)=#; polytomy2=[0,3] #=(1,2,(3,4))=#; polytomy3=[0,4] #=((1,2),3,4)=#
    polytomy4=[1,3] #=(1,(2,3,4))=#; polytomy5=[2,4] #=((1,2,3),4)=#; polytomy6=[1,2] #=(1,(2,3),4)=#

    for node in q.mrca
        ij=node[1]; ik=node[2]; il=node[3]; jk=node[4]; jl=node[5]; kl=node[6]
        if ij==ik==il==jk==jl==kl push!(q.symtype,rand(polytomy1))#@debug "At least one of the quartets have hard polytomy for all taxa."
        elseif ij==ik==il==jk==jl#=!==kl=# push!(q.symtype,rand(polytomy2))#@debug "At least one of the quartets have hard polytomy for first two taxa."
        elseif ik==il==jk==jl==kl#=!==ij=# push!(q.symtype,rand(polytomy3))#@debug "At least one of the quartets have hard polytomy for last two taxa."
        elseif ij==ik==il==jl==kl#=!==jk=# push!(q.symtype,rand(polytomy6))#@debug "At least one of the quartets have hard polytomy for center two taxa."
        elseif ik==il==jk==jl push!(q.symtype,0)#Symmetric Tree
        elseif ij==ik==il
            if jk==jl==kl push!(q.symtype,rand(polytomy4))#@debug "At least one of the quartets have hard polytomy for last three taxa."
            elseif jl==kl push!(q.symtype,1)#Asymmetric Type I
            else push!(q.symtype,3)#Asymmetric Type III
            end
        elseif ij==ik
            if jl==kl
                if jk==ij push!(q.symtype,rand(polytomy5))#@debug "At least one of the quartets have hard polytomy for first three taxa."
                else push!(q.symtype,2)#Asymmetric Type II
                end
            end
        else push!(q.symtype,4) #=Asymmetric Type IV=# 
        end
    end
    return q
end


function uniqtaus(q::nquartets)
    for n in q.mrca
        if q.symtype[1]==0 t1=tauNum(n[6]); t2=tauNum(n[1]); t3=tauNum(n[3]) #((1,2),(3,4)) [x,r,r,r,r,y]
        elseif q.symtype[1]==1 t1=tauNum(n[4]); t2=tauNum(n[5]); t3=tauNum(n[3]) #(1,((2,3),4)) [r,r,r,x,y,y]
        elseif q.symtype[1]==2 t1=tauNum(n[4]); t2=tauNum(n[1]); t3=tauNum(n[3]) #((1,(2,3)),4) [x,x,r,y,r,r]  
        elseif q.symtype[1]==3 t1=tauNum(n[6]); t2=tauNum(n[4]); t3=tauNum(n[3]) #(1,(2,(3,4))) [r,r,r,x,x,y]
        elseif q.symtype[1]==4 t1=tauNum(n[1]); t2=tauNum(n[4]); t3=tauNum(n[3]) #(((1,2,),3),4) [x,y,r,y,r,r]
        end
        push!(q.ntau,[t1,t2,t3])
    end
end

function DictionaryP(p::Phylip)
    d=Dict(p.nametaxa[n]=>p.counttaxa[n] for n=1:length(p.nametaxa))
        return d
end

function DictionaryQ(q::Nquartets)
    d=Dict(q.leafname[n]=>q.leafnumber[n] for n=1:length(q.leafname))
        return d
end

function DictionaryN(q::Nquartets,p::Phylip)
    dp=DictionaryP(p)
    dq=DictionaryQ(q)

    if length(dp)!==length(dq)
        println("The number of individuals in Phylip and the tree does not match.")
    end

    dn=Dict{Int64,Int64}(dq[x]=>dp[x] for x in q.leafname)
    
    return dn
end

function pquartet(q::Nquartets,dn::Dict)
    for x in q.nquartet
        #get quartet numbers equivalent in the phylip file
        i=dn[x.quartet[1][1]]; j=dn[x.quartet[1][2]]; k=dn[x.quartet[1][3]]; l=dn[x.quartet[1][4]]
        push!(x.tquartet,[i,j,k,l])
    end
end

function binaryIndexforQuartet(q::Nquartets)
    for x in q.nquartet
        i=bitstring(Int8(x.tquartet[1][1])); j=bitstring(Int8(x.tquartet[1][2]))
        k=bitstring(Int8(x.tquartet[1][3])); l=bitstring(Int8(x.tquartet[1][4]))
        push!(x.indexNET,[i*j*k*l])
    end
end

function moveSPcounts(q::Nquartets,p::Phylip)
    dn=DictionaryN(q,p)
    pquartet(q,dn)
    binaryIndexforQuartet(q)
    for x in q.nquartet
        for m in 1:length(p.allquartet)
            if x.indexNET[1]==p.index[m] push!(x.mspcountsNET,p.spcounts[m]) end
        end 
    end
end

#=function findGamma(subtree::HybridNetwork,net::HybridNetwork)
    gamma=Float64[]
    if net.numHybrids==0
        gamma=[1.0]
    else
        for n in net.edge
            if n.hybrid==true
                x=GetParent(n).number
                for t in subtree.edge
                    z=GetParent(t).number
                    if x==z
                        push!(gamma,n.gamma)
                        break
                    else 
                        continue
                    end
                end
            else
                continue
            end
        end
    end
    return prod(gamma)
end=#

"""
    printQuartets

gegejaga
"""
printQuartets(x) = printQuartets(stdout::IO, x)
function printQuartets(io::IO, qq::Array)
    count=0
    for q in qq
        count+=1
        println(io, "Gamma for parental tree $count: $(q.gamma)")
        println(io, "#\t q [i,j,k,l]\tmrca [ij,ik,il,jk,jl,kl]\ttau#\t\ttype")
        for q in q.nquartet 
            println("$(q.number)\t$(q.quartet[1])\t$(q.mrca[1])\t$(q.ntau[1])\t$(q.symtype)") 
        end
        println("\n")
    end
end