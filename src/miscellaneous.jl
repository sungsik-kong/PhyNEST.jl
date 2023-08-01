#Written by Sungsik Kong 2021-2022
#Last updated by Sungsik Kong 2023
"""
    greet()

Displays a greeting with citation information. No argument needed. 

## Example
```@jldoctest
julia> greet()
Thank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.
Please report bugs or make suggestions to https://groups.google.com/g/phynest-users.
If you conduct an analysis using PhyNEST, please cite:
Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood.
Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.
```
"""
function greet()
    #now=Dates.now()
    println("""
    Thank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.
    Please report bugs or make suggestions to https://groups.google.com/g/phynest-users.
    If you conduct an analysis using PhyNEST, please cite:
    Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood.
    Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.
    """)
end


"""
    get_average(i::Array)

Compute average of the values in an array. Did not want to add Statistics dependency for this...
"""
function get_average(i::Array)
    average=sum(i)/length(i)
    return average
end



"""
    checkDEBUG()

Just a simple function to check if debugging message is being displayed. If debugging message is not being displayed, try:
- ENV["JULIA_DEBUG"] = "PhyNEST" or 
- ENV["JULIA_DEBUG"] = "All".
To turn it off, try:
- ENV["JULIA_DEBUG"] = " ".
"""
function checkDEBUG()
    println("Function is being executed. Do you see a message on next line? If not try ?checkDEBUG().")
    @debug "Debugging message is being displayed."
end


















































"""
    correct_outgroup(net::HybridNetwork, outgroup::AbstractString)

    Check if the outgroup is indeed the outgroup of the network
"""
function correct_outgroup(net::HybridNetwork, outgroup::AbstractString)
    rooted_with_outgroup=false
    root_node_number=net.root
    root=net.node[root_node_number]
    if length(root.edge)==2
        attached_edge_1_to_root=root.edge[1]
        attached_edge_2_to_root=root.edge[2]

        child1=GetChild(attached_edge_1_to_root)
        child2=GetChild(attached_edge_2_to_root)

        if child1.name==outgroup
            rooted_with_outgroup=true
            return rooted_with_outgroup
        elseif child2.name==outgroup
            rooted_with_outgroup=true
            return rooted_with_outgroup
        else
            rooted_with_outgroup=false
            return rooted_with_outgroup
        end
    end
    return rooted_with_outgroup
end




function LRT(p::Phylip)
    tree=readTopology("(1,(2,(3,4)));")
    net=readTopology("(1,((2,(3)#H6:::0.5),(4,#H6:::0.5)));")

    stat0,top0=do_optimization(tree,p)
    stat1,top2=do_optimization(net,p)

    likelihood0=stat0.minimum
    likelihood1=stat1.minimum

    lrt=-2 * (likelihood0-likelihood1)

    

    return lrt 

end


























Dstat(outgroup::String,p::Phylip)=Dstat(outgroup::String,p::Phylip,0.05,false)
Dstat(outgroup::String,p::Phylip,pval::Float64)=Dstat(outgroup::String,p::Phylip,pval::Float64,false)
function Dstat(outgroup::String,p::Phylip,pval::Float64,displayall::Bool)
    
    outgroup1=0
    ndist=Normal(0,1)
    quartet=[]
    sitepattern=[]

    for n in 1:length(p.nametaxa)
        
        if p.nametaxa[n]==outgroup
            outgroup1=p.counttaxa[n]
            break
        end
    end

    outgroup1!==0 || error("Cannot find the specified outgroup '$outgroup' in Phylip input file.")

    for n in 1:length(p.allquartet)
        if p.allquartet[n][1]==outgroup1
            push!(quartet,p.allquartet[n])
            push!(sitepattern,[p.spcounts[n][7],p.spcounts[n][9]])
        else
            continue
        end
    end

        for n in 1:length(sitepattern)
            ABAB=sitepattern[n][1]
            ABBA=sitepattern[n][2]
            d=(ABBA-ABAB)/(ABBA+ABAB)
            z=d/(2*sqrt((0.25/(ABBA+ABAB))))
            println([ABAB, ABBA, d, z])
            
            pv=1-cdf(ndist,z)
            if displayall==true
                println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv")
            else
                if pv<=(pval)
                    println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv <= SIGNIFICANT!")
                else
                    println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv")
                end
            end

        end

    end


###Need to deal with avg_obs=num_obs/(1*1*1*1)#(self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
function HyDe(outgroup::String,p::Phylip; p_value=0.05::Float64, filter=true::Bool)
    outgroup1=0
    quartet=[]
    sitepattern=[]
    num_obs = p.seqleng
     
    

    for n in 1:length(p.nametaxa)
        if p.nametaxa[n]==outgroup
            outgroup1=p.counttaxa[n]
            break
        end
    end

    outgroup1!==0 || error("Cannot find the specified outgroup '$outgroup' in Phylip input file.")

    for n in 1:length(p.allquartet)
        if p.allquartet[n][1]==outgroup1
            push!(quartet,p.allquartet[n])
            push!(sitepattern,[p.spcounts[n][4],p.spcounts[n][7],p.spcounts[n][9],p.spcounts[n][3],p.spcounts[n][2],p.spcounts[n][6]])
        else
            continue
        end
    end


    for n in 1:length(sitepattern)

        p9 = (sitepattern[n][3]+0.05)/num_obs#ABBA
        p7 = (sitepattern[n][2]+0.05)/num_obs#ABAB
        p4 = (sitepattern[n][1]+0.05)/num_obs#AABB

        avg_obs=num_obs/(1*1*1*1)#(self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
        avobs=avg_obs
        
        #z-score
        obs_invp1 = avobs * (p9 - p7)
        obs_invp2 = avobs * (p4 - p7)
        obs_var_invp1 = avobs * p9 * (1 - p9) + avobs * p7 * (1 - p7) + 2 * avobs * p9 * p7
        obs_var_invp2 = avobs * p4 * (1 - p4) + avobs * p7 * (1 - p7) + 2 * avobs * p4 * p7
        obs_cov_invp1_invp2 = -1 * avobs * p9 * p4 + avobs * p9 * p7 + avobs * p7 * p4 + avobs * p7 * (1 - p7)
        ratio = obs_invp2 / obs_invp1;
        GH_ts = ((obs_invp1) * (ratio) / sqrt(obs_var_invp1 * (ratio^2) - 2.0 * obs_cov_invp1_invp2 * ratio + obs_var_invp2))

        temp = -99999.9
        if p7 > p9 && p7 < p4
            GH_ts=temp
        elseif GH_ts > -99999.9 && GH_ts < 99999.9
            GH_ts=GH_ts
        else
            GH_ts=temp
        end
        my_z=GH_ts


        #gamma
        _c_num   = avg_obs * (p9 - p7)
        _c_denom = avg_obs * (p4 - p7)
        _c       = _c_num / _c_denom
        gamma = _c / (1 + _c)


        #p-val
        a1 =  0.254829592;
        a2 = -0.284496736;
        a3 =  1.421413741;
        a4 = -1.453152027;
        a5 =  1.061405429;
        px  =  0.3275911;
        sign = 1;
        if my_z < 0 sign = -1 end
        z = abs(my_z) / sqrt(2.0);
        t = 1.0 / (1.0 + px * z);
        y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-z * z);

        p_val=1.0 - (0.5 * (1.0 + sign * y))


        newQ=[p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]
        if filter
            if p_val<(p_value/3)
                println("$newQ, gamma: $gamma, z-score: $GH_ts, p-val: $p_val")
            end
        else
            println("$newQ, gamma: $gamma, z-score: $GH_ts, p-val: $p_val")
        end
    end

end



















# add a leaf 
using PhyloNetworks
#q=readTopology("(1,2,(3,4));")
function printEverything(t::HybridNetwork) PhyloNetworks.printEverything(t) end


function add_a_leaf(top::HybridNetwork, new_leaf_name::AbstractString)
    treeset=HybridNetwork[]
    #Create a branch (u,v) where 
    #u=leaf with leaf number=n+1; where n=number of leaves
    #v=smallest node number-1
    #println("The new leaf name is set as '$new_leaf_name'")
    n=length(top.leaf)
    #println("There are $n leaves in the topology..")
    u=n+1
    smallest_internal_node_number=-2
    for nods in top.node
        if nods.number < smallest_internal_node_number
            smallest_internal_node_number=nods.number
        end
    end
    v=smallest_internal_node_number-1

    #make the node objects for u and v
    node_u=PhyloNetworks.Node(u,true)
    node_u.name=new_leaf_name
    node_v=PhyloNetworks.Node(v,false)
    
    #make an edge uv that links u and v
    #make the edge number as current edges + 1 because we are adding an edge
    edge_uv=PhyloNetworks.Edge(top.numEdges+1)
    edge_uv.node=[node_u,node_v]
    edge_uv.length=-1.0
    #add information to nodes u and v that they are linked to edge_uv
    push!(node_u.edge,edge_uv)
    push!(node_v.edge,edge_uv)
    #push!(t.node,node_u)
    #push!(t.node,node_v)
    #println("=====node_u, node_v, edge_uv=====")
    #println(node_u)
    #println(node_v)
    #println(edge_uv)
    
    #for a selected edge (u',v')
    #   add (u,v) - push 2 nodes
    #   create 2 edges (u',v), (v,v')
    #   remove (u',v')

    for i in 1:length(top.edge)
        t=deepcopy(top)
        e=t.edge[i]
        #println("=====tree before modify, selected edge=====")
        #printEverything(t)
        #println(e)

        #creating two edges that disects the selected edge
        numEdges=length(t.edge)
        e1=deepcopy(e)
        e1.number=numEdges+1
        e1.node[2]=node_v #tail of one of the edges is v
        e2=deepcopy(e)
        e2.number=numEdges+2        
        e2.node[1]=node_v #head of one of the edges is v  
        #println("=====selected edge, modified edge 1, modified edge 2=====")
        #println(e, e1, e2)      

        #so far we need to remove e and push e1 and e2
        #also we need to add node_u, node_v, and edge_uv
        deleteat!(t.edge, i)
        push!(t.edge,e1)
        push!(t.edge,e2)
        push!(t.node,node_u)
        push!(t.node,node_v)
        push!(t.edge,edge_uv)
        #println("=====tree after modify=====")
        #printEverything(t)

        #now everything that is needed is out There
        #so reoragnize the edge numbers connected to each node 
        #since some are deleted and some are newly added
        for vertex in t.node
            relevantedges=PhyloNetworks.Edge[]
            thenumber=vertex.number
            for branch in t.edge
                parent=PhyloNetworks.getParent(branch)
                parent=parent.number
                child=PhyloNetworks.getChild(branch)
                child=child.number
                if parent == thenumber
                    push!(relevantedges,branch)
                elseif child == thenumber
                    push!(relevantedges,branch)
                else
                    continue
                end
            end
            vertex.edge=relevantedges
        end    
        #println("=====tree alsmot final=====")
        #printEverything(t)
        t0=readTopology(PhyloNetworks.writeTopologyLevel1(t))
        push!(treeset,t0)
    end
    
    return treeset
end


function add_n(tree::HybridNetwork, new_leaf_name::AbstractString)
    treeset=HybridNetwork[]
    push!(treeset,tree)
    col=add_n(treeset, new_leaf_name)
    return col
end

function add_n(trees::Array, new_leaf_name::AbstractString)
    treeset=HybridNetwork[]
    for eachtree in trees
        set=add_a_leaf(eachtree, new_leaf_name)
        for each in set
            push!(treeset,each)
        end
    end

    return treeset
end

#make a function that pipelines all - return object or newick
function add_one(trees::Array; new_leaf_name=""::AbstractString)
    treeset=HybridNetwork[]
    for eachtree in trees
        set=add_n(eachtree, new_leaf_name)
        for each in set
            push!(treeset,each)
        end
    end
    return treeset
end
function add_one(top::HybridNetwork; new_leaf_name=""::AbstractString)
    if isempty(new_leaf_name) new_leaf_name=length(top.leaf)+1 end
    new_leaf_name=(string(new_leaf_name))
    collection=add_n(top,new_leaf_name)
    return collection
end

#given a topology, switch leaf labels and return the list
function flip_leaves(tree::HybridNetwork)
    set=HybridNetwork[]
    #Replace species names?
    return set
end