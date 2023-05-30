#Written by Sungsik Kong 2021-2022
#Last updated by Sungsik Kong 2023

###---create quartet object only with topology---###
"""
    get_quartets(net::HybridNetwork; round_gamma_digits=5::Int64)
    get_quartets(net::HybridNetwork, p::Phylip; round_gamma_digits=5::Int64, default_theta=0.0001::Float64)

Extracts information for every quartet extractable from the topology and store it as a Network object. \\
When the Phylip object is also provided, we fill more slots including the theta, branch lengths for \\
eac quartet using the moment estimator, average branch lengths using the estimated branch lengths.

### Input
- `net`                   PhyloNetwork.HybridNetwork object for a topology [mandatory]\\
- `p`                     Phylip object\\
- `round_gamma_digits`    Number of digits to show for the inhertiance probability in each quartet\\
- `default_theta          Theta is being estimated when both topology and data are given, however, \\
                        when the estimated theta is beyond the expected range (i.e., 0.00001 < theta < 0.1)\\
                        we use user-specified theta (default=0.0001) to estimate branch lengths using \\
                        the moment estimator.
"""
function get_quartets(net::HybridNetwork; round_gamma_digits=5::Int64)
    N=Network()
    N.leafname=get_leaf_name(net)
    N.leafnumber=get_leaf_number(net)
    N.gamma=get_gamma(net)
    numind=length(N.leafnumber)
    quartet_number=1
    displayed_trees=displayedTrees(net,0.0) #displayedTrees function from PhyloNetworks
    numtrees=length(displayed_trees)
    displayed_tree_number=1
        
    while displayed_tree_number <= numtrees
        quartet_number_in_this_tree=1
        current_tree=displayed_trees[displayed_tree_number]
        current_tree_leafnumber=get_leaf_number(current_tree)
        all_quartets_in_this_disaplyed_tree=list_all_quartets(current_tree_leafnumber)

        while quartet_number_in_this_tree <= (binomial(numind,4))
            this_quartet=quartets(quartet_number,displayed_tree_number)
            this_quartet.gamma=round(N.gamma[displayed_tree_number],digits=round_gamma_digits)
            this_quartet.quartet=all_quartets_in_this_disaplyed_tree[quartet_number_in_this_tree]
            this_quartet.mrca=get_most_recent_common_ancestors(displayed_trees[displayed_tree_number],this_quartet)
            this_quartet.symtype=identify_sym_type(this_quartet)
            this_quartet.ntau=get_unique_tau_labels(this_quartet)
            push!(N.quartet,this_quartet)    

            quartet_number_in_this_tree+=1
            quartet_number+=1 
        end
        displayed_tree_number+=1
    end
    return N
end

###---create quartet object with topology AND phylip information (i.e. fill more information)---###
#using theta=0.0001 by default
function get_quartets(net::HybridNetwork, p::Phylip; 
                        round_gamma_digits=5::Int64, default_theta=0.0001::Float64)
    N=get_quartets(net;round_gamma_digits=round_gamma_digits)
    for each_quartet in N.quartet
        dictionary_Combined=dictionary_combined(N,p)
        each_quartet.tquartet=get_matching_quartet(each_quartet,dictionary_Combined)
        each_quartet.mspcountsNET=transfer_site_pattern_frequencies(each_quartet,p)
    end
    estimated_theta=default_theta
    estimated_theta=get_start_theta(N)
    if 0.00001 < estimated_theta < 0.1
        N.theta=estimated_theta
    else
        N.theta=estimated_theta 
        #@warn("Expected 0.00001 < estimated_theta < 0.1, but received $estimated_theta. Now N.theta=$default_theta.")
    end

    for each_quartet in N.quartet
        each_quartet.momestlength=momentEstimate(each_quartet,estimated_theta)
    end
    average_momest=get_average_moment_branch_length(N)
    for each_quartet in N.quartet
        each_quartet.average_mom_est_bl=(average_momest[each_quartet.ntau[1]],average_momest[each_quartet.ntau[2]],average_momest[each_quartet.ntau[3]])
    end
    return N
end

"""
    get_leaf_name(net::HybridNetwork)

Extracts name of each leaf in the topology from the HybridNetwork object.
"""
function get_leaf_name(net::HybridNetwork)
    leafname=String[]
    for tip in net.leaf
        (tip.leaf) && push!(leafname, tip.name)
    end
    length(unique(leafname))==length(leafname) || error("Some taxa appears more than once in the input network.")    
    return leafname
end

"""
    get_leaf_number(net::HybridNetwork)

Extracts number of each leaf in the topology from the HybridNetwork object. \\
This number is different from the sequence number stored in Phylip object.
"""
function get_leaf_number(net::HybridNetwork)
    leafnumber=Int64[]
    for tip in net.leaf
        (tip.leaf) && push!(leafnumber, tip.number)
    end
    return leafnumber
end

"""
    list_all_quartets(leaf_number::Array)

Using the leaf numbers given in HybridNetworks, obtains all quartet extract-able from the topology.
"""
function list_all_quartets(leaf_number::Array)
    all_quartets=[]
    numind=length(leaf_number) 
    numind > 0 || error("There is no leaf specified in the input network")

    for i in (1:numind)
        for j in (i+1:numind)
            for k in (j+1:numind)
                for l in (k+1:numind)
                    push!(all_quartets,(leaf_number[i],leaf_number[j],leaf_number[k],leaf_number[l]))
                end
            end
        end
    end

    return all_quartets
end

"""
    GetChild(edge::Edge)

Pasted function from PhyloNetworks. See documentation for PhyloNetworks.GetChild
"""
function GetChild(edge::Edge) edge.node[edge.isChild1 ? 1 : 2] end

"""
    GetParent(edge::Edge)

Pasted function from PhyloNetworks. See documentation for PhyloNetworks.GetParent
"""
function GetParent(edge::Edge) edge.node[edge.isChild1 ? 2 : 1] end

"""
    childnode(parentnode::Int64,network::HybridNetwork)

Returns child node of the given the parent node number stored in HybridNetwork. In case there are two children (i.e., tree node), it spits out `one of the two` nodes.
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
    parentnode(childnode::Int64,network::HybridNetwork)

Returns parent node of the given the child node number stored in HybridNetwork. In case there are two parents (i.e., hybrid node), it spits out `one of the two` nodes.
"""
function parentnode(childnode::Int64,network::HybridNetwork)
    for edge in network.edge
        if childnode==GetChild(edge).number
            parentnode=GetParent(edge).number
        return parentnode
        end
    end
end

"""
    mrca(Node1::Int64,Node2::Int64,net::HybridNetwork)

Returns the most recent common ancestor node for two given nodes. By nodes, I mean the node number stored in HybridNetwork.
"""
function mrca(Node1::Int64,Node2::Int64,net::HybridNetwork)
    root=net.node[net.root].number
    root!==nothing || @error("Root cannot be identified in the input topology")
    route1=Int64[]
    route2=Int64[]
    push!(route1,Node1)
    push!(route2,Node2)

    while Node1!==root 
        Node1=parentnode(Node1,net)
        Node1!==nothing || @error("No parent exists for $Node1 in the input topology.")
        typeof(Node1)==Int64 || @error("Parent node of the node $Node1 is weird.")
        push!(route1,Node1)
    end

    while Node2!==root
        Node2=parentnode(Node2,net)
        Node2!==nothing || @error("No parent exists for $Node2 in the input topology.")
        typeof(Node2)==Int64 || @error("Parent node of the node $Node2 is weird.")
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

"""
    get_most_recent_common_ancestors(net::HybridNetwork,each_quartet::quartets)    

Identifies the most recent common ancestor for every pair of leaves and \\
arrange those into an array that contains six to identify the tree shape of the quartet in identify_sym_type().
"""
function get_most_recent_common_ancestors(net::HybridNetwork,each_quartet::quartets)    
    each_quartet_member=each_quartet.quartet
    i=each_quartet_member[1] 
    j=each_quartet_member[2]
    k=each_quartet_member[3] 
    l=each_quartet_member[4]

    ij=mrca(i,j,net)
    ik=mrca(i,k,net)
    il=mrca(i,l,net)
    jk=mrca(j,k,net)
    jl=mrca(j,l,net)
    kl=mrca(k,l,net)

    return (ij,ik,il,jk,jl,kl)
end

"""
    identify_sym_type(each_quartet::quartets)

Looking at the six most recent common ancestors for every pair of leaves obtained from get_most_recent_common_ancestors() \\
and recorded in quartets.mrca, this function identifies if the quartet is symmetric (tyep 0) or asymmetric; and if asymmetric\\
it further identifies which one of the four asymmetric quartet (type 1, 2, 3, 4) it is. Although unpreferred, if a quartet turns out to be\\
a polytomy, we randomly select any one of the five possible topologies. For example, any one of types 0, 1, 2, 3, or 4 will be selected \\
if the quartet turns out to be a hard polytomy (i.e., (i,j,k,l);) where the six mrcas will be something like [x,x,x,x,x,x].
"""
function identify_sym_type(each_quartet::quartets)
    polytomy1=[0,1,2,3,4] #=(1,2,3,4)=#; polytomy2=[0,3] #=(1,2,(3,4))=#; polytomy3=[0,4] #=((1,2),3,4)=#
    polytomy4=[1,3] #=(1,(2,3,4))=#; polytomy5=[2,4] #=((1,2,3),4)=#; polytomy6=[1,2] #=(1,(2,3),4)=#

    mrca_node=each_quartet.mrca
    ij=mrca_node[1]; ik=mrca_node[2]; il=mrca_node[3]; jk=mrca_node[4]; jl=mrca_node[5]; kl=mrca_node[6]
    if ij==ik==il==jk==jl==kl sym_type=rand(polytomy1)#@debug "At least one of the quartets have hard polytomy for all taxa."
    elseif ij==ik==il==jk==jl#=!==kl=# sym_type=rand(polytomy2)#@debug "At least one of the quartets have hard polytomy for first two taxa."
    elseif ik==il==jk==jl==kl#=!==ij=# sym_type=rand(polytomy3)#@debug "At least one of the quartets have hard polytomy for last two taxa."
    elseif ij==ik==il==jl==kl#=!==jk=# sym_type=rand(polytomy6)#@debug "At least one of the quartets have hard polytomy for center two taxa."
    elseif ik==il==jk==jl sym_type=0#Symmetric Tree
    elseif ij==ik==il
        if jk==jl==kl sym_type=rand(polytomy4)#@debug "At least one of the quartets have hard polytomy for last three taxa."
        elseif jl==kl sym_type=1#Asymmetric Type I
        else sym_type=3#Asymmetric Type III
        end
    elseif ij==ik
        if jl==kl
            if jk==ij sym_type=rand(polytomy5)#@debug "At least one of the quartets have hard polytomy for first three taxa."
            else sym_type=2#Asymmetric Type II
            end
        end
    else sym_type=4 #=Asymmetric Type IV=# 
    end

    return sym_type
end

"""
    tauNum(x::Int64)

Convert the node number assigned in PhyloNetworks to the format used in PhyNEST.
"""
function tauNum(x::Int64) return abs(x+1) end

"""
    backtauNum(x::Int64)

Back-convert the node number assigned in PhyNEST to the format used in PhyloNetworks.
"""
function backtauNum(x::Int64) return -1*(x+1) end

"""
    get_unique_tau_labels(each_quartet::quartets)

For the array that contains six mrac node numbers, it selects three of those that are unique from each other\\
and returns them in the order of increasing node age.
"""
function get_unique_tau_labels(each_quartet::quartets)
    q=each_quartet; n=q.mrca
    if q.symtype[1]==0 t1=tauNum(n[6]); t2=tauNum(n[1]); t3=tauNum(n[3]) #((1,2),(3,4)) [x,r,r,r,r,y]
    elseif q.symtype[1]==1 t1=tauNum(n[4]); t2=tauNum(n[5]); t3=tauNum(n[3]) #(1,((2,3),4)) [r,r,r,x,y,y]
    elseif q.symtype[1]==2 t1=tauNum(n[4]); t2=tauNum(n[1]); t3=tauNum(n[3]) #((1,(2,3)),4) [x,x,r,y,r,r]  
    elseif q.symtype[1]==3 t1=tauNum(n[6]); t2=tauNum(n[4]); t3=tauNum(n[3]) #(1,(2,(3,4))) [r,r,r,x,x,y]
    elseif q.symtype[1]==4 t1=tauNum(n[1]); t2=tauNum(n[4]); t3=tauNum(n[3]) #(((1,2,),3),4) [x,y,r,y,r,r]
    end
    return (t1,t2,t3)
end

"""
    dictionary_phylip(p::Phylip)

Creates a dictionary that collects the name of taxa provided in the phylip format to the taxa number stored in Phylip object.
"""
function dictionary_phylip(p::Phylip)
    dictionary_Phylip=Dict(p.nametaxa[n]=>p.counttaxa[n] for n=1:length(p.nametaxa))
        return dictionary_Phylip
end

"""
    dictionary_topology(N::Network)

Creates a dictionary that collects the names of leaf provided in the topology to the leaf number stored in Network object \\
(that are extracted from HybridNetwork.PhyloNetworks).
"""
function dictionary_topology(N::Network)
    dictionary_Topology=Dict(N.leafname[n]=>N.leafnumber[n] for n=1:length(N.leafname))
        return dictionary_Topology
end

"""
    dictionary_combined(N::Network,p::Phylip)

Creates a dictionary that collects the taxa number in Phylip and leaf number in Network by combinging the two dictinoaries\\
produced in dictionary_phylip() and dictionary_topology(). 
"""
function dictionary_combined(N::Network,p::Phylip)
    dictionary_Phylip=dictionary_phylip(p)
    dictionary_Topology=dictionary_topology(N)

    length(dictionary_Phylip)==length(dictionary_Topology) || @error("The number of taxa in Phylip does not match with the number of leaves in tooplogy.")

    dictionary_combined=Dict{Int64,Int64}(dictionary_Topology[x]=>dictionary_Phylip[x] for x in N.leafname)
    
    return dictionary_combined
end

"""
    function get_matching_quartet(each_quartet::quartets,dictionary_combined::Dict)

Using the dictionary created in dictionary_combined(), it converts the quartet number that was given during parsing the topoogy\\
to the matching number that is given and stored in the Phylip object. Depending on the topology, this step may not matter.
"""
function get_matching_quartet(each_quartet::quartets,dictionary_combined::Dict)
    #get quartet numbers equivalent in the phylip file
    i=dictionary_combined[each_quartet.quartet[1]]; 
    j=dictionary_combined[each_quartet.quartet[2]]; 
    k=dictionary_combined[each_quartet.quartet[3]]; 
    l=dictionary_combined[each_quartet.quartet[4]]
    return (i,j,k,l)
end

"""
    function transfer_site_pattern_frequencies(each_quartet::quartets,p::Phylip)

Using the quartet number that is obtained in get_matching_quartet(), it transfers the matching site pattern frequencies stored\\
in the Phylip object. This will be the `relevant' observed data that is used for our analysis for that particular topology.
"""
function transfer_site_pattern_frequencies(each_quartet::quartets,p::Phylip)
    for quartet_number in 1:length(p.allquartet)
        phylip_quartet_tuple=(p.allquartet[quartet_number][1],p.allquartet[quartet_number][2],p.allquartet[quartet_number][3],p.allquartet[quartet_number][4]) 
        if each_quartet.tquartet==phylip_quartet_tuple
            return p.spcounts[quartet_number]
        end
    end
end

"""
    gamArray

Returns the inheritance probability for each parental tree of the network.
"""
function get_gamma(net::HybridNetwork)
    given_gammas=Float64[]
    parental_trees_with_gamma=displayedTrees(net,0.0,multgammas=true)
    for each_tree in parental_trees_with_gamma
        tree_gamma=1.0
        each_edge=each_tree.edge
        for edge in each_edge
            tree_gamma*=edge.gamma
        end
            push!(given_gammas,tree_gamma)
    end
    return given_gammas
end

"""
    printQuartets(N::Network)

Prints relevant information of the quartets extracted from the topology.
"""
###---Display quartet information---###
printQuartets(N) = printQuartets(stdout::IO, N)
function printQuartets(N::Network)
    number_of_displayed_trees=length(N.gamma)
    all_quartets=N.quartet
    #basic information
    println("Summary of the quartets extracted from the topology")
    println("Number of displayed trees: $(number_of_displayed_trees)")
    println("Quartet\tDisplayed tree\tQuartet\t\tGamma\tTau numbers\tQuartet type")
    for q in all_quartets
        println("$(q.number)\t$(q.displayed_tree)\t\t$(q.quartet)\t$(q.gamma)\t$(q.ntau)\t$(q.symtype)")
    end
end
