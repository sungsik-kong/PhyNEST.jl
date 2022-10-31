#written by Sungsik Kong 2021-2022

# 6.Network Search
const Move2Int = Dict{Symbol,Int}(:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5,:nni=>6)
const Nature = Dict{Int,String}(0=>"rootNode",1=>"internal node",2=>"leaf",3=>"hybrid node")

"""
Updating topology Parameters

nature
updateTopolEdges
updateTopolGammas
updateTopolParams
"""

"""
    nature(Node::Int64,net::HybridNetwork)

Function to determine the nature (root, internal tree node, leaf, hybrid node) of node
given the node number and Hybridnetwork.

### Input
- **`Node`**       An integer node number\\
- **`net`**        A tree/network in HybridNetwork\\
"""
function nature(Node::Int64,net::HybridNetwork)
    #nature 0=root;1=tree node;2=leaf;3=hybrid node
    nat=-1
    numNodes=length(net.node)
    rootNodeNum=net.root
    root=net.node[rootNodeNum].number
    for node in 1:numNodes
        if net.node[node].number==Node
            if net.node[node].leaf nat=2
            elseif net.node[node].hybrid nat=3
            elseif Node==root nat=0
            else nat=1
            end
        end
    end
    nat>=0 || @debug "Cannot determine the nature of node $Node."
    #@debug "Node $(Node) is a $(Nature[nat])"
    return nat
end

"""
    updateTopolEdges(net::HybridNetwork,taus::Array)

Function to update estimates of speciation times to Topology. Sometimes the estimates can be weird, like containing
NaN (Not a Number), and this may complicate things. Therefore, we warn the users about when this happens. 

### Input
- **`net`**        A tree/network in HybridNetwork\\
- **`taus`**       An array of estimates of speciation times\\
"""
function updateTopolEdges(net::HybridNetwork,taus::Array)
    !isnan(sum(taus)) || @debug "One or more estimates of speciation times may not be a number. Estimates are: $taus.)"
    #we try to prevent the possiblilty of above error in [updateTopolParams] before executing this function.
    for e in net.edge
        parentnode=GetParent(e).number
        childnode=GetChild(e).number
        natp=nature(parentnode,net)
        natc=nature(childnode,net)
        #weird natp and natc are prevented in [nature]
        if natp!==3 && natc!==3 #if both parent and child are not hybrid nodes
            if childnode>0 #if childnode is a leaf
                natc==2 || @debug "Childnode is expected to be a leaf, instead we got $(Nature[natc]). "
                e.length=taus[tauNum(parentnode)]
            else #or else the childnode must be a internal tree node
                natc==1 || @debug "Childnode is expected to be an internal tree node, instead we got $(Nature[natc]). "
                e.length=taus[tauNum(parentnode)]-taus[tauNum(childnode)]
            end
        end
    end
    @debug "updateTopolEdges: Estimated τ $taus successfully updated to topology"
    return net
end

"""
    updateTopolGammas(net::HybridNetwork,gammas::Array,retedge::Array)

Function to update estimates of inheritance probabilities to Topology. Sometimes the estimates can be weird, like containing
NaN (Not a Number), and this may complicate things. Therefore, we warn the users about when this happens. 

### Input
- **`net`**        A tree/network in HybridNetwork\\
- **`gammas`**     An array of estimates of inheritance probabilities\\
- **`retedge`**    An array of (u,v) for reticulation edge in topology that has gamma.
- For example, in the case of retedge=[[-5,5],[-8,2]], gammas=[[0.5],[0.1]], the edge (-5,5) has gamma of 0.5 and the edge (-8,2) has gamma of 0.1.\\
"""
function updateTopolGammas(net::HybridNetwork,gammas::Array,retedge::Array)
    !isnan(sum(gammas)) || @debug "One or more estimates of inheritance probability may not be a number. Estimates are: $gammas.)"
    for e in net.edge
        if (e.hybrid)
            parentnode=GetParent(e).number
            childnode=GetChild(e).number
            for i in 1:length(retedge)
                if childnode==retedge[i][2]
                    if parentnode==retedge[i][1]
                        e.gamma=gammas[i]
                    else
                        e.gamma=(1-gammas[i])
                    end
                end
            end
        end
    end
    @debug "Estimated ɣ $gammas on edge u,v=$retedge successfully updated to topology"
    return net
end

"""
    updateTopolParams(net::HybridNetwork,taus::Array,gammas::Array,retedge::Array)

Function to use estimates of parameters from optimization and update onto topology if the numbers seem to be
useful (i.e., none of the estimates is NaN). Otherwise, it does not update 
"""
function updateTopolParams(net::HybridNetwork,taus::Array,gammas::Array,retedge::Array)
    if !isnan(sum(taus)) && !isnan(sum(gammas))
        @debug "Both estimates of τ and ɣ seems to be useful. Updating parameters onto topology."
        net=updateTopolEdges(net,taus)
        net=updateTopolGammas(net,gammas,retedge)
    elseif !isnan(sum(taus)) 
        net=updateTopolEdges(net,taus)
        @debug "Some estimates of ɣ is not a number (NaN). Not updating ɣ on topology."
    elseif !isnan(sum(gammas))
        net=updateTopolGammas(net,gammas,retedge)
        @debug "Some estimates of τ is not a number (NaN). Not updating τ on topology."
    else 
        @debug "Estimates of τ and ɣ seems to be not useful. Not updating parameters onto topology."
    end
    return net
end

#removes branch lengths from a topology
"""
    fixNegatBranch(net::HybridNetwork)

Function to remove all existing branch lengths by setting it to -1.0. This is to prevent a network to have negative or NaN branch length
that will cause error during optimization.

### Input
- **`net`**        A tree/network in HybridNetwork\\
"""
function fixNegatBranch(net::HybridNetwork)
    suc=false
    for e in net.edge 
        suc=false
        e.length=-1.0 
        suc=true
    end
    if (suc) @debug "Successfully removed all τ on topology." end
    return net
end

#set all gamma to 0.5
"""
    fixNegatGamma(net::HybridNetwork)

Function to initialize all existing gamma by setting it to 0.5. This is to prevent a network to have negative or NaN or reticulation
that does not add up to 1.0 which will cause error during optimization.

### Input
- **`net`**        A tree/network in HybridNetwork\\
"""
function fixNegatGamma(net::HybridNetwork)
    suc=false
    for e in net.edge
        e.gamma >= 0.0 && e.gamma<=1.0 || @debug "ɣ must be between 0 and 1 but got $(e.gamma) for edge $(e.number) ."
        while e.gamma !== 1.0 
            suc=false
            e.gamma=0.5 
            suc=true
            break end
    end
    if (suc) @debug "Successfully set ɣ to 0.5 for $(net.numHybrids) reticulations." end
    return net
end

"""
    fixWierdos(net::HybridNetwork)

This function takes a network in HybridNetwork and removes all branch lengths and set all gamma to 0.5

### Input
- **`net`**        A tree/network in HybridNetwork\\
"""
function fixWierdos(net::HybridNetwork)
    net=fixNegatGamma(net)
    net=fixNegatBranch(net)
end

"""
    multStartT(net::HybridNetwork,outgroup::AbstractString)

Function to make an NNI modification on the topology so the search can be initiated in different starting 
points. This could be very useful if we have low confidence in the starting topology.
In brief, using the input network, it unroots the topology, conduct NNI modification, removes branch lengths 
and set gamma to 0.5, and root the topology using the ougroup provided. We assume this should be (almost always) successful.

### Input
- **`net`**        A tree/network in HybridNetwork\\
- **`outgroup`**   Outgroup taxa in a string\\
"""
function multStartT(net::HybridNetwork,outgroup::AbstractString)
    suc=false
    unrootN=readTopologyUpdate(writeTopologyLevel1(net))#unroot the net
    while (!suc)
        suc=NNIRepeat!(unrootN,10)
    end
    nStartT=fixWierdos(unrootN) #unrooted, check if all gamma is (0,1); if all branch legnth is >0
    @suppress begin nStartT=readTopology(writeTopologyLevel1(nStartT,outgroup)) end #Roots the network
    if (suc) @debug "Successfully modified the starting topology using NNI." 
    else @debug "Modying the starting topology using NNI unsuccessful." end
    return nStartT
end

"""
    function hillClimb(startT::HybridNetwork,p::Phylip,hmax::Integer,outgroup::String,maxcount::Integer,
    nfail::Integer,NumIter::Integer,paramprint::Bool,logfile::IO,display::Bool)

    Executes Hill-Climbing Algorithm.
"""
function hillClimb(startT::HybridNetwork,p::Phylip,hmax::Integer,outgroup::String,maxcount::Integer,
        nfail::Integer,NumIter::Integer,paramprint::Bool,logfile::IO,display::Bool)
    #set some necessary stuff for search
    steps = 0
    failures = 0
    movescount = zeros(Int,18)#1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
    movesfail = zeros(Int,6)#count of failed moves for current topology
    Nmov = zeros(Int,6) 
    stillmoves = true
    #get mpl for the starting network
    @suppress begin currT=readTopology(writeTopologyLevel1(startT)) end
    fixWierdos(currT)
    currTLogLik,taus,gammas,theta2,retedge=Optimization(currT,p,NumIter,paramprint)
    currT=updateTopolParams(currT,taus,gammas,retedge)
    #set newT=currT
    newT=deepcopy(currT)
    #begin search
    while(steps<maxcount && failures < nfail && stillmoves)
        steps+=1 #tracking number of topologies evaluated 
        calculateNmov!(newT,Nmov)# function to calculate Nmov, number max of tries per move2int [0,0,0,0,0,0]=>[122,25,25,4,10000,42]
        move = whichMove(newT,hmax,movesfail,Nmov) #selects which move to take
        if move != :none
            accepted=false
            newT=@suppress begin readTopologyUpdate(writeTopologyLevel1(newT,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,newT,true,steps,10, movescount,movesfail,false) #unrooted, make modification on newT accroding to move
            newT=fixWierdos(newT) #unrooted, check if all gamma is (0,1); if all branch legnth is >0
            newT=@suppress begin readTopology(writeTopologyLevel1(newT,outgroup)) end#Roots the network
            newTLogLik,taus,gammas,theta2,retedge=Optimization(newT,p,NumIter,false)
            newT=updateTopolParams(newT,taus,gammas,retedge)
            if(newTLogLik<currTLogLik)
                accepted=true
            else
                accepted=false
            end

            if(accepted)
                currT=newT
                currTLogLik=newTLogLik
                failures = 0
                movescount[Move2Int[move]+12] += 1
                movesfail = zeros(Int,6)#count of failed moves for current topology back to zero 
            else
                newT=currT
                failures += 1
                movesfail[Move2Int[move]] += 1
            end
        else
            stillmoves=false
        end
    end

    #determine the reason for termination
    termination=0
    if steps==maxcount termination=1
    elseif failures==nfail termination=2
    end
    #print reasons for termination
    str = 
"The search terminated at step $steps and at $(failures)th consecutive failures.
Summary of each move:
Insertion of reticulation edge: $(movescount[1]) proposed, $(movescount[13]) accepted.
Tail move of reticulation edge: $(movescount[2]) proposed, $(movescount[14]) accepted. 
Head move of reticulation edge: $(movescount[3]) proposed, $(movescount[15]) accepted.
Change the direction of reticulation edge: $(movescount[4]) proposed, $(movescount[16]) accepted.
Deletion of reticulation edge: $(movescount[5]) proposed, $(movescount[17]) accepted.
Nearest-neighbor interchange (NNI): $(movescount[6]) proposed, $(movescount[18]) accepted.
On the current topology, $(sum(movesfail)) moves were made, including $(sum(movesfail)-failures) unsuccessful moves.\n"
    if termination==1 
        str *= 
"Terminated because it reached the maximum number of steps (current maxcount=$maxcount). 
It is recommended to increase the number of maxcount and rerun the analysis.\n"
    elseif termination==2 
        str *= 
"Terminated because it reached the maximum number of failures (current nfail=$nfail).\n"
    else 
        str *= 
"Terminated although it neither reached the maximum number of steps nor failures,
possibly because there was no more move to make.\n"
    end
    if(display) print(str) end
    write(logfile,str); flush(logfile)   

    return currTLogLik,currT

end

function burnin!(startT::HybridNetwork, p::Phylip, hmax::Integer, outgroup::String, burninn::Int64,NumIter::Integer)

    BurnIn=Float64[]

    steps = 0
    movescount = zeros(Int,18)
    movesfail = zeros(Int,6)
    Nmov = zeros(Int,6)
    stillmoves = true

    currT=startT
    mpl,taus,gammas,theta2,retedge=Optimization(currT,p,NumIter,false)
    currTLogLik=mpl

    newT=deepcopy(currT)

    while(length(BurnIn)<burninn)
        steps+=1

        calculateNmov!(newT,Nmov)# function to calculate Nmov, number max of tries per move2int
        move = whichMove(newT,hmax,movesfail,Nmov)

        if move != :none
            newT=readTopologyUpdate(writeTopologyLevel1(newT)) #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,newT,true,steps,10, movescount,movesfail,false) #unrooted, make modification on newT accroding to move
            newT=fixWierdos(newT) #unrooted, check if all gamma is (0,1); if all branch legnth is >0
            @suppress begin newT=readTopology(writeTopologyLevel1(newT,outgroup)) end #Roots the network
            newTLogLik,taus,gammas,theta2,retedge=Optimization(newT,p,NumIter,false)

            delta=abs(currTLogLik-newTLogLik)

            currTLogLik=newTLogLik
            if !isnan(delta) && !iszero(delta) && !isinf(delta)
                push!(BurnIn,delta)
            end

        else
            stillmoves=false
        end
    end

    return BurnIn
end

#burnin!(x,p,1,"5",10,1000)

function maxburnin(startT::HybridNetwork, p::Phylip, hmax::Integer, outgroup::String, burninn::Int64,NumIter::Integer)
    BurnIn=burnin!(startT,p,hmax,outgroup,burninn,NumIter)
    if isempty(BurnIn)
        BurnIn=burnin!(startT,p,hmax,outgroup,burninn,NumIter)
    else       
        x=findmax(BurnIn)[1] 
        return x
    end
end

#maxburnin(x,p,1,"5",10,1000)

function nfails(p::Phylip, prob::Float64) 
    n=p.numtaxa
    x=(log(prob)/log(1-(1/(n-2)))/n)
    x=round(Int64,x)
    return x
end

#i think a=c=0.5 kind of works better...

function simAnneal!(startT::HybridNetwork, p::Phylip, hmax::Integer, outgroup::String, burninn::Int64,maxcount::Integer,
                    k::Integer,ProbabilityOfNotSelectingNode::Float64,cons::Float64,alph::Float64,NumIter::Integer,
                    paramprint::Bool,logfile::IO,display::Bool)
    
    #stating some necessaries...
    ktrees=[]
    #kktrees=[]
    #determining the number of failured based on the prob.
    nfail=nfails(p,ProbabilityOfNotSelectingNode)
    steps = 0
    failures = 0
    movescount = zeros(Int,18)
    movesfail = zeros(Int,6)
    Nmov = zeros(Int,6)
    stillmoves = true
    ci=0.0

    #dealing with the starting tree
    currT=deepcopy(startT)
    mpl,taus,gammas,theta2,retedge=Optimization(currT,p,NumIter,paramprint)
    currT=updateTopolParams(currT,taus,gammas,retedge)
    currTLogLik=mpl
    newT=currT
    cT=writeTopologyLevel1(currT)
    push!(ktrees,(currTLogLik,cT))
    #push!(kktrees,(currTLogLik,cT))

    #burn-in 
    str="Running $burninn runs of burn-in..." 
    if(display) print(str) end
    write(logfile,str); flush(logfile)
    u=maxburnin(currT, p, hmax, outgroup, burninn, NumIter) 
    l=currTLogLik
    b=cons/(((1-alph)*p.numtaxa)+((alph)*(log.(l)/p.seqleng)))
    str="Complete 
(Cooling schedule)U=$u
(Cooling schedule)beta=$b\n"
    if(display) print(str) end
    write(logfile,str); flush(logfile)
    
    while(steps < maxcount && failures < nfail && stillmoves)
        #1
        steps+=1
        i=steps-1
        calculateNmov!(newT,Nmov)# function to calculate Nmov, number max of tries per move2int
        move = whichMove(newT,hmax,movesfail,Nmov)

        if move != :none
            accepted=false
            #newT=readTopologyUpdate(writeTopologyLevel1(newT))
            newT=@suppress begin readTopologyUpdate(writeTopologyLevel1(newT,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,newT,true,steps,10, movescount,movesfail,false) 
            newT=fixWierdos(newT)
            @suppress begin newT=readTopology(writeTopologyLevel1(newT,outgroup)) end
            
            #2
            mpl,taus,gammas,theta2,retedge=Optimization(newT,p,NumIter,paramprint)
            newT=updateTopolParams(newT,taus,gammas,retedge)
            newTLogLik=mpl
            if newTLogLik<=currTLogLik && !isnan(newTLogLik) && !isnan(currTLogLik)
                accepted=true
            elseif newTLogLik>currTLogLik && !isnan(newTLogLik) && !isnan(currTLogLik)
                ci=u/(1+(i*b))
                prob=exp((-1*(newTLogLik-currTLogLik))/ci)
                
                #println("L1-L0=$(newTLogLik-currTLogLik)")
                #println("ci=$ci")
                #println("t=$prob\n")

                items=[true,false]
                weigh=[prob,1-prob]
                accepted=sample(items, Weights(weigh))
            else
                accepted=false
            end
            if(accepted)
                
                nT=writeTopologyLevel1(newT)
                ktrees=sort!(ktrees, by = x -> x[1])
                flag=false
                for e in ktrees
                    if (newTLogLik in e[1])==false
                        flag=false 
                    else
                        flag=true 
                        break
                    end
                end
                if flag==true
                    failures += 1
                    movesfail[Move2Int[move]] += 1
                else
                    failures = 0
                    movesfail = zeros(Int,6) 
                    movescount[Move2Int[move]+12] += 1
                end
                
                if length(ktrees)<k
                    if steps==1
                        push!(ktrees,(currTLogLik,cT))
                        #push!(kktrees,(currTLogLik,currT))
                    else
                        if flag==false
                            push!(ktrees,(newTLogLik,nT))
                            #push!(kktrees,(newTLogLik,newT))
                        end
                    end
                elseif length(ktrees)==k && flag==false
                    #@suppress begin nT=writeTopologyLevel1(newT,outgroup) end
                    if newTLogLik<ktrees[k][1]
                        ktrees[k]=(newTLogLik,nT)
                    end
                end
                currT = newT
                currTLogLik=newTLogLik
            else
                newT = currT
                failures += 1
                movesfail[Move2Int[move]] += 1
            end
        else
            stillmoves=false;
        end
    end
    #kktrees=sort!(kktrees, by = x -> x[1])
    ktrees=sort!(ktrees, by = x -> x[1])

#    println(kktrees)#[mpl, HybridNetwork Objects]
    #println(ktrees)#[mpl, newick]

    rank=0
    str =
"Rank   Pseudolikelihood    Network\n"
    for i in ktrees
        rank+=1
    str *=
"$rank\t$(round(i[1],digits=5))         $(i[2])\n"
    end
    
    if(rank<k) str *="The total length of best trees can be shorter than k.\n" end
    str *="Speciation times for some newicks may not have updated if estimates are weird (e.g., NaN).\n"
    if(display) print(str) end
    write(logfile,str); flush(logfile)   

    mpl=ktrees[1][1]
    bestnet=ktrees[1][2]
     #determine the reason for termination
     termination=0
     if steps==maxcount termination=1
     elseif failures==nfail termination=2
     end
     #print reasons for termination
     str = 
"The search terminated at step $steps and at $(failures)th consecutive failures and (Cooling schedule)ci=$ci.
Summary of each move:
Insertion of reticulation edge: $(movescount[1]) proposed, $(movescount[13]) accepted.
Tail move of reticulation edge: $(movescount[2]) proposed, $(movescount[14]) accepted. 
Head move of reticulation edge: $(movescount[3]) proposed, $(movescount[15]) accepted.
Change the direction of reticulation edge: $(movescount[4]) proposed, $(movescount[16]) accepted.
Deletion of reticulation edge: $(movescount[5]) proposed, $(movescount[17]) accepted.
Nearest-neighbor interchange (NNI): $(movescount[6]) proposed, $(movescount[18]) accepted.
On the current topology, $(sum(movesfail)) moves were made, including $(sum(movesfail)-failures) unsuccessful moves.\n"
    if termination==1 
        str *= 
"Terminated because it reached the maximum number of steps (current maxcount=$maxcount). 
It is recommended to increase the number of maxcount and rerun the analysis.\n"
    elseif termination==2 
        str *= 
"Terminated because it reached the maximum number of failures (current nfail=$nfail).\n"
    else 
        str *= 
"Terminated although it neither reached the maximum number of steps nor failures,
possibly because there was no more move to make.\n"
    end
    if(display) print(str) end
    write(logfile,str); flush(logfile)   

    return mpl,bestnet

end


"""
    PhyNE!(startT::HybridNetwork,inputFile::Phylip,outgroup::String)

Estimates the network or tree to fit observed site pattern frequencies stored in a `Phylip` object,
using composite likelihood. A level-1 network is assumed. The search begins from topolgoy `startT`,
which can be a tree or a network (with less than `hmax` reticulation nodes). Usually, `startT` is
expected to be a tree estiamted from data, but can also be randomly generated. The topology is rooted
using the `outgroup` provided. This must be identical to the outgroup sequence provided in the 
sequence alignment. By default, `PhyNE!` will search the network space using simulated annealing assuming 
hmax=1.

There are many optional arguments (values in parenthesis are defaults):
- `hillclimbing(=false)` Select hill climbing search
- `hmax(=1)` Maximum number of reticulations in the final network
- `nruns(=10)` Number of independent runs
- `nniStartT=false` 

- `maxcount=1000` Number of steps to terminate the search
- `nfail=75` Number of consecutive failures to terminate the search

Simulated annealing stuff:
- `burninn=25`
- `k=10`
- `ProbabilityOfNotSelectingNode=9.5e-45`
- `cons=0.9`
- `alph=0.8`

Optimization stuff:
- `NumIter=1000`
- `paramprint=false`

Miscellaneous:
- `timestamp=false`
- `filename=""`
- `display=false`s
"""
function PhyNE!(startT::HybridNetwork,inputFile::Phylip,outgroup::String; 
    hmax=1::Integer, nruns=5::Integer,maxcount=100000::Integer,NumIter=1000::Integer,
    hillclimbing=false::Bool, nfail=75::Integer, burninn=25::Int64, k=10::Integer,
    ProbabilityOfNotSelectingNode=9.5e-45::Float64,cons=0.9::Float64,alph=0.8::Float64,
    paramprint=false::Bool,timestamp=false::Bool,nniStartT=false::Bool,
    filename=""::AbstractString,display=false::Bool)

    #current time
    t=Dates.now()
    times=Dates.format(t,"mmddHHMM")

    #filename
    if isempty(filename)
        filename="PhyNE"
        if(hillclimbing) filename=filename*(".hc")
        else filename=filename*(".sa") end
    end
    #timestamp for filename
    if(timestamp) filename=filename*(".")*(string(times)) end

    #restate some stuff
    p=inputFile
    T=deepcopy(startT)
    i=0 #count nruns
    errr=false
    net=[]

    #Logging some information in prior to actual Analysis
    if(hillclimbing) algo=String("Hill-climbing") else algo=String("Simulated Annealing") end
    log=string(filename,".log")
    logfile=open(log,"w")
    str =
"PhyNE: Estimating Maximum Pseudolikelihood Phylogenetic Network
╔═══╦╗─────╔═╗─╔╗
║╔═╗║║─────║║╚╗║║
║╚═╝║╚═╦╗─╔╣╔╗╚╝╠══╗
║╔══╣╔╗║║─║║║╚╗║║║═╣
║║──║║║║╚═╝║║─║║║║═╣
╚╝──╚╝╚╩═╗╔╩╝─╚═╩══╝
───────╔═╝║ 
───────╚══╝   
Analysis start: $(Dates.format(t, "yyyy-mm-dd at HH:MM:SS"))
Input filename: $(p.filename)
Number of sequences: $(p.numtaxa) 
Sequence length: $(p.seqleng)
Starting Topology: $(writeTopologyLevel1(startT))
Outgroup specified for rooting: $outgroup
Number of maximum reticulation(s): $hmax
The maximum number of iterations for each optimization: $NumIter
Search algorithm selected: $algo
The maximum number of steps during search: $maxcount\n"
if(hillclimbing==false) str*="Alpha: $alph; Cons: $cons\n" end
str*="\nInitiating $nruns iterations...\n"
    if(display) print(str) end
    write(logfile,str)
    flush(logfile)

    #begin analysis
    if(hillclimbing)
        while i < nruns
            i+=1
            if(nniStartT) T=multStartT(T,outgroup) end
            str=
"($i/$nruns) Searching for the best network using the hill-climbing algorithm...
Starting topology modified to $(writeTopologyLevel1(T))\n"
            if(display) print(str) end
            write(logfile,str)
            flush(logfile)
            try 
                mpl,BestNet=hillClimb(T,p,hmax,outgroup,maxcount,nfail,NumIter,paramprint,logfile,display)
                bT=writeTopologyLevel1(BestNet)
                push!(net,[mpl,BestNet])
                str = 
"The best network found in this run: $bT
-Log Pseudolikelihood: $mpl \n\n"
                if(display) print(str) end
                write(logfile,str); flush(logfile)    
            catch(error)
                errr=true
                print("error")
            end
        end
    else
        while i < nruns
            i+=1
            if(nniStartT) T=multStartT(T,outgroup) end
            str=
"($i/$nruns) Searching for the best network using the simulated annealing algorithm...
Starting topology modified to $(writeTopologyLevel1(T))\n"
            if(display) print(str) end
            write(logfile,str)
            flush(logfile)
            try
                mpl,BestNet=simAnneal!(startT,p,hmax,outgroup,burninn,maxcount,k,ProbabilityOfNotSelectingNode,cons,alph,NumIter,paramprint,logfile,display)
                #bT=writeTopologyLevel1(BestNet)
                push!(net,[mpl,BestNet])
                str = 
"The best network found in this run: $BestNet
-Log Pseudolikelihood: $mpl \n\n"
                if(display) print(str) end
                write(logfile,str); flush(logfile)    
            catch(error)
                errr=true
                print("error")
            end
        end
    end

    #writing outfile
    out=string(filename,".out")
    outfile=open(out,"w")

    net=sort!(net, by = x -> x[1])
    e=Dates.now()#Dates.format(now(), "yyyy-mm-dd at HH:MM:SS")
    deltatsec= (e-t) / Millisecond(1) * (1 / 1000)
    if(hillclimbing)
        bT=writeTopologyLevel1(net[1][2])
        dscopeT=writeTopologyLevel1(net[1][2],di=true)
    else
        net1=readTopology(net[1][2])
        bT=writeTopologyLevel1(net1)
        dscopeT=writeTopologyLevel1(net1,di=true)
    end

    str = "Summary:\n"; write(logfile,str); flush(logfile)      
    if(hillclimbing) str = "The best network found from $nruns runs using the hill-climbing algorithm\n"; 
    else str = "The best network found from $nruns runs using the simulated annealing algorithm\n"; end
    write(outfile,str); flush(outfile)      
    str = 
    "MPL network: \n$(bT)
    Dendroscope: \n$(dscopeT)
    -Log Pseudolikelihood: $(net[1][1]).\nend\n"
    write(logfile,str); flush(logfile)      
    write(outfile,str); flush(outfile)
    str =
    "―――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    MPL network: \n$(bT)
    Analysis complete: $(Dates.format(e, "yyyy-mm-dd at HH:MM:SS"))
    Time Taken: $(round(deltatsec,digits=5)) second(s)\nend\n\n"; print(stdout,str)
    write(logfile,str); flush(logfile)      

    return net[1][2]
end
