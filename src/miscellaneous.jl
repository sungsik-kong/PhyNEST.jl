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
    Dstat(outgroup::String, p::Phylip;
        pval=0.05::Float64, 
        displayall=false::Bool)

Conducts Patterson's D-statistic test. The result prints site pattern frequencies ABAB and ABBA used to compute 
the D-statistic, Z-score, and the p-value for each quartet tested. Significance is marked with an asterisk.
Function `showall(df)` can be subsequently used to show all rows.

## Mandatory arguments
- `outgroup`     Name of the outgroup taxa
- `p`   The `Phylip` object

## Optional arguments
- `pval       (default=0.05)` Alpha level for significance
- `display_all (default=false)` If set as `true`, the function print test results for every quartet. By default, it only prints those quartets where signficance was found.

## Example
```@jldoctest
julia> p=readPhylip("sample_n4h1.txt")
julia> df=Dstat("4",p)
Tip: if neccessary, use showall(df) function to see all the rows.
2×10 DataFrame
 Row │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat     Zscore   pvalue   significance
     │ String    String  String  String  Int64  Int64  Float64   Float64  Float64  String
─────┼──────────────────────────────────────────────────────────────────────────────────────────
   1 │ 4         3       2       1        1427   7852  0.692424  66.6995      0.0  *
   2 │ 4         1       2       3        1427   7836  0.691892  66.5908      0.0  *

julia> df=Dstat("4",p,display_all=true)
Tip: if neccessary, use showall(df) function to see all the rows.
6×10 DataFrame
Row │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat        Zscore      pvalue    significance
    │ String    String  String  String  Int64  Int64  Float64      Float64     Float64   String
────┼─────────────────────────────────────────────────────────────────────────────────────────────────
  1 │ 4         3       1       2        7852   1427  -0.692424    -66.6995    1.0
  2 │ 4         3       2       1        1427   7852   0.692424     66.6995    0.0       *
  3 │ 4         1       3       2        7836   1427  -0.691892    -66.5908    1.0
  4 │ 4         1       2       3        1427   7836   0.691892     66.5908    0.0       *
  5 │ 4         2       3       1        7836   7852   0.00101989    0.127743  0.449176
  6 │ 4         2       1       3        7852   7836  -0.00101989   -0.127743  0.550824
```
"""
function Dstat(outgroup::String, p::Phylip; pval=0.05::Float64, display_all=false::Bool)

    dict=dictionary_phylip(p)
    ndist=Normal(0,1)
    res=[]

    #Reassure the provided outgroup indeed exists in the data
    outgroup_presence=issubset([outgroup],p.nametaxa)
    (outgroup_presence) || error("Outgroup $outgroup does not exist in data.")
    ourgroup_id=dict["$outgroup"]

    #ABBA-BABA test
    for n in 1:length(p.allquartet)
        if p.allquartet[n][1]==ourgroup_id 
            out=p.allquartet[n][1]
            t1=p.allquartet[n][2]
            t2=p.allquartet[n][3]
            t3=p.allquartet[n][4]

            ABAB=p.spcounts[n][7]
            ABBA=p.spcounts[n][9]

            d=(ABBA-ABAB)/(ABBA+ABAB)
            z=d/(2*sqrt((0.25/(ABBA+ABAB))))
            pv=1-cdf(ndist,z)

            if pv<=(pval) ast="*" else ast="" end

            if (display_all)
                push!(res,[p.nametaxa[out],p.nametaxa[t1],p.nametaxa[t2],p.nametaxa[t3],ABAB, ABBA, d, z, pv, ast])
            else
                if ast=="*"
                    push!(res,[p.nametaxa[out],p.nametaxa[t1],p.nametaxa[t2],p.nametaxa[t3],ABAB, ABBA, d, z, pv, ast])
                end
            end
        else
            continue
        end
    end

    #prepare to print results in a cleaner way using DataFrame
    df=DataFrame(outgroup=String[],
                taxa1=String[],
                taxa2=String[],
                taxa3=String[],
                ABAB=Int[],
                ABBA=Int[],
                Dstat=Float64[],
                Zscore=Float64[],
                Pvalue=Float64[],
                sig=String[])
    for result in res
        push!(df, result)
    end            
   
    println("Tip: if neccessary, use showall(df) function to see all the rows.")

    return df
end

"""
    showall(df::DataFrame)    

Print all rows of the DataFrame object.    
"""
function showall(df::DataFrame) CSV.show(df,allrows=true)   end
























function Dstatsearchquartets()
end

"""

"""
###Need to deal with avg_obs=num_obs/(1*1*1*1)#(self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
###i.e., multiple individuals / sp

#check if 1 ind res matches
#check if multi ind yields the same spfreq
#check res in various sim
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











































"###---Adding a leaf on a topology: Misc. function designed for a task with Tang at UWM---###"
# add a leaf 
#q=readTopology("(1,2,(3,4));") #unrooted quartet

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


#transforms a topology into a adjacent matrix
#a topology is transformed to n by n matrix where n=number of nodes
#0 if there is no linkage between the nodes; 1 if there is a linkage
function top2adjmat4tip(top::HybridNetwork)
    number_of_nodes=length(top.node)
    dict=Dict()

    dict["1"]=1
    dict["2"]=2
    dict["3"]=3
    dict["4"]=4
    dict["5"]=5
    dict["6"]=6
    
    length(dict)==number_of_nodes || error("Number of nodes and length of dictionary does not match.")
    
    adjmat=zeros(Int64, number_of_nodes, number_of_nodes)
    for vertex in top.node
        edges=vertex.edge
        for link in edges
            x=dict[vertex.name]
            y1=dict[GetParent(link).name]
            y2=dict[GetChild(link).name]
            adjmat[x,x]=1
            adjmat[x,y1]=1
            adjmat[x,y2]=1
        end
    end
    #println(adjmat)
    return adjmat
end

function top2adjmat5tip(top::HybridNetwork)
    number_of_nodes=length(top.node)
    dict=Dict()

    dict["1"]=1
    dict["2"]=2
    dict["3"]=3
    dict["4"]=4
    dict["5"]=5
    dict["6"]=6
    dict["7"]=7
    dict["8"]=8
    
    length(dict)==number_of_nodes || error("Number of nodes and length of dictionary does not match.")
    
    adjmat=zeros(Int64, number_of_nodes, number_of_nodes)
    for vertex in top.node
        edges=vertex.edge
        for link in edges
            x=dict[vertex.name]
            y1=dict[GetParent(link).name]
            y2=dict[GetChild(link).name]
            adjmat[x,x]=1
            adjmat[x,y1]=1
            adjmat[x,y2]=1
        end
    end
    #println(adjmat)
    return adjmat
end