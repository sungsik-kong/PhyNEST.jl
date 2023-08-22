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
    Dstat(p::Phylip, outgroup::String;
        alpha=0.05::Float64, 
        displayall=false::Bool,
        writecsv=false::Bool, 
        filename=""::AbstractString)

Conducts Patterson's D-statistic test. The result prints site pattern frequencies ABAB and ABBA used to compute 
the D-statistic, Z-score, and the p-value for each quartet tested. Significance is marked with an asterisk.
Function `showall(df)` can be subsequently used to show all rows.

## Mandatory arguments
- `p`   The `Phylip` object
- `outgroup`     Name of the outgroup taxa

## Optional arguments
- `alpha       (default=0.05)` Alpha level for significance
- `display_all (default=false)` If set as `true`, the function print test results for every quartet. By default, it only prints those quartets where the signficance was found.
- `writecsv (default=false)` If `true`, the result is stored in `.csv` file in the working directory
- `filename` Specifies `.csv` file name if `writecsv=true`. If unspecified, the result is stored as `Dstat-out.csv`

## Example
```@jldoctest
julia> p=readPhylip("sample_n4h1.txt")
julia> df=Dstat(p,"4")
Tip: if neccessary, use showall(df) function to see all the rows.
2×10 DataFrame
 Row │ outgroup  taxa1   taxa2   taxa3   ABAB   ABBA   Dstat     Zscore   pvalue   significance
     │ String    String  String  String  Int64  Int64  Float64   Float64  Float64  String
─────┼──────────────────────────────────────────────────────────────────────────────────────────
   1 │ 4         3       2       1        1427   7852  0.692424  66.6995      0.0  *
   2 │ 4         1       2       3        1427   7836  0.691892  66.5908      0.0  *

julia> df=Dstat(p,"4",display_all=true)
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
function Dstat(p::Phylip, outgroup::String; 
                alpha=0.05::Float64, 
                display_all=false::Bool, 
                writecsv=false::Bool, 
                filename=""::AbstractString)

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

            if pv<=(alpha) ast="*" else ast="" end

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
                significance=String[])
    for result in res
        push!(df, result)
    end            
   
    #write csv
    if filename=="" filename="Dstat-out" end
    if (writecsv)
        CSV.write("$filename.csv",df) 
        println("The results are stored as $filename.csv in the working directory.")    
    end

    println("Tip: if neccessary, use function showallDF(df) to see all the rows.")

    return df
end

"""
    showallDF(df::DataFrame)    

Print all rows of the DataFrame object using the package `CSV`.    
"""
function showallDF(df::DataFrame) CSV.show(df,allrows=true)   end




"""
    HyDe(p::Phylip, outgroup::AbstractString; 
        alpha=0.05::Float64, 
        display_all=true::Bool, #filter
        map=""::AbstractString,
        writecsv=false::Bool, 
        filename=""::AbstractString


Conducts HyDe: Hybrid Detectioo using Phylogenetic invariants. See Blischak et al., (2018) (https://doi.org/10.1093/sysbio/syy023) and the manual (https://hybridization-detection.readthedocs.io) for more information. 
This function replicates `run_hyde.py` in the original python package. The map file, a two-column table with individual names in the first column and the name of the population that it belongs to in the second column
(see example below) is not necessary to start the analysis. If the map file is not provided, each sequnece in the data is assumed to represent distinct species. 

Map file is a simple text file where each line contains the name of the sequencea and the assignment of the sequence to a species, delimited by a tab. See an example below.

## Mandatory arguments
- `p`   The `Phylip` object        
- `outgroup`     Name of the outgroup taxa. Even when there are muliple outgroup individuals, but if their assignment to the species is provided in the `map` file, simply specify one of the outgroup species as listed in the alignment.


## Optional arguments
- `alpha       (default=0.05)` Alpha level for significance
- `display_all (default=false)` If set as `true`, results are also filtered based on whether there is significant evidence for hybridization.
- `map (default=no map file)`   Specify a map file, if available. 
- `writecsv (default=false)` If `true`, the result is stored in `.csv` file in the working directory
- `filename` Specifies `.csv` file name if `writecsv=true`. If unspecified, the result is stored as `HyDe-out.csv`

## Example
```@jldoctest
julia> HyDe(p,"5",display_all=true)
Tip: if neccessary, use function showallDF(df) to see all the rows.
24×11 DataFrame
 Row │ outgroup  P1      Hybrid  P2      AABB   ABAB   ABBA   Gamma         Zscore         Pvalue    significance
     │ String    String  String  String  Int64  Int64  Int64  Float64       Float64        Float64   String
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 5         4       3       1       18168   2573   2611    0.00243076       0.528406  0.298609
   2 │ 5         4       3       2       23742   2022   2064    0.00192997       0.657666  0.255376
   3 │ 5         4       1       2       23703   2069   2069    0.0         -99999.9       1.0
   4 │ 5         3       1       2        8005   8057   1991    0.9915          -0.412067  0.659855
   5 │ 5         4       1       3       18168   2611   2573   -0.00244861  -99999.9       1.0
   6 │ 5         3       4       1        2573  18168   2611    0.49939       -216.327     1.0
   7 │ 5         3       1       4        2573   2611  18168    1.00245         -0.527118  0.700944
   8 │ 5         1       4       3        2611  18168   2573    0.50061       -216.327     1.0
   9 │ 5         1       3       4        2611   2573  18168    0.997569         0.528406  0.298609
  10 │ 5         4       2       3       23742   2064   2022   -0.00194121  -99999.9       1.0
  11 │ 5         3       4       2        2022  23742   2064    0.499516      -339.45      1.0
  12 │ 5         3       2       4        2022   2064  23742    1.00194         -0.656395  0.744215
  13 │ 5         2       4       3        2064  23742   2022    0.500484      -339.45      1.0
  14 │ 5         2       3       4        2064   2022  23742    0.99807          0.657666  0.255376
  15 │ 5         4       2       1       23703   2069   2069    0.0         -99999.9       1.0
  16 │ 5         1       4       2        2069  23703   2069    0.5           -336.307     1.0
  17 │ 5         1       2       4        2069   2069  23703  NaN                0.0       0.5
  18 │ 5         2       4       1        2069  23703   2069    0.5           -336.307     1.0
  19 │ 5         2       1       4        2069   2069  23703  NaN                0.0       0.5
  20 │ 5         3       2       1        8005   1991   8057    0.502152        47.6571    0.0       *
  21 │ 5         1       3       2        8057   8005   1991    1.00872     -99999.9       1.0
  22 │ 5         1       2       3        8057   1991   8005    0.497848        47.6571    0.0       *
  23 │ 5         2       3       1        1991   8005   8057   -0.00872191      -0.408534  0.658559
  24 │ 5         2       1       3        1991   8057   8005    0.00849951      -0.412067  0.659855

julia> HyDe(p,"5",display_all=false)
Tip: if neccessary, use function showallDF(df) to see all the rows.
2×11 DataFrame
 Row │ outgroup  P1      Hybrid  P2      AABB   ABAB   ABBA   Gamma     Zscore   Pvalue   significance
     │ String    String  String  String  Int64  Int64  Int64  Float64   Float64  Float64  String
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 5         3       2       1        8005   1991   8057  0.502152  47.6571      0.0  *
   2 │ 5         1       2       3        8057   1991   8005  0.497848  47.6571      0.0  *

julia> HyDe(p,"5",map="map.txt",display_all=false)
Map file [map.txt] provided.
Tip: if neccessary, use function showallDF(df) to see all the rows.
2×11 DataFrame
 Row │ outgroup  P1      Hybrid  P2      AABB   ABAB   ABBA   Gamma     Zscore   Pvalue   significance
     │ String    String  String  String  Int64  Int64  Int64  Float64   Float64  Float64  String
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ sp5out    sp3     sp2     sp1     15841   3418  15909  0.501365  49.4337      0.0  *
   2 │ sp5out    sp1     sp2     sp3     15909   3418  15841  0.498635  49.4337      0.0  *
```

where the `map.txt` used in above example is:
```
shell> cat map.txt
5	sp5out
4	sp5out
3	sp3
1	sp1
2	sp2
```
"""
function HyDe(p::Phylip, outgroup::AbstractString; 
                alpha=0.05::Float64, 
                display_all=false::Bool, #filter
                map=""::AbstractString,
                writecsv=false::Bool, 
                filename=""::AbstractString)
    
    res=[]
    data=[]
    outgroups=[]
    num_obs = p.seqleng
    havemap=!(isempty(map))
    
    #Reassure the provided outgroup indeed exists in the data
    outgroup_presence=issubset([outgroup],p.nametaxa)
    (outgroup_presence) || error("Outgroup $outgroup does not exist in data.")
    
    #if mapfile provided reorganize data[] for the analysis    
    #presence/absence of the map file
    if havemap 
        println("Map file [$map] provided.")
        mapd=readdlm(map, '\t', AbstractString, '\n')
        map_nrows=size(mapd,1)
        map_nrows==p.numtaxa || error("Number of taxa in the alignment does not match with the map file.")
        mapdict=Dict{AbstractString,AbstractString}()
            for i in 1:map_nrows
                mapdict[mapd[i,1]]=mapd[i,2]
            end
            outgroup_id=mapdict["$outgroup"]    
            for i in 1:map_nrows
                if mapd[i,2]==outgroup_id
                    push!(outgroups,mapd[i,1])
                else continue
                end
            end
    else
        push!(outgroups,outgroup)
    end

    #get the quartets and relevant site patterns for the analysis
    for n in 1:length(p.allquartet)
        for outgroup_ids in outgroups
            if p.nametaxa[p.allquartet[n][1]]==outgroup_ids
                if havemap
                    t1=mapdict[p.nametaxa[p.allquartet[n][1]]]
                    t2=mapdict[p.nametaxa[p.allquartet[n][2]]]
                    t3=mapdict[p.nametaxa[p.allquartet[n][3]]]
                    t4=mapdict[p.nametaxa[p.allquartet[n][4]]]
                else
                    t1=p.nametaxa[p.allquartet[n][1]]
                    t2=p.nametaxa[p.allquartet[n][2]]
                    t3=p.nametaxa[p.allquartet[n][3]]
                    t4=p.nametaxa[p.allquartet[n][4]]                 
                end
                
                if t1==t2 continue
                elseif t1==t3 continue
                elseif t1==t4 continue
                else 
                    push!(data,[p.allquartet[n][1],
                                p.allquartet[n][2],
                                p.allquartet[n][3],
                                p.allquartet[n][4],
                                p.spcounts[n][4],
                                p.spcounts[n][7],
                                p.spcounts[n][9],
                                t1,
                                t2,
                                t3,
                                t4])
                end
            end
        end
    end
    
    #summarize and add up according to the map file
    if !(isempty(map))
        newdata=[]
        for i in 1:length(data)-1
            for j in i+1:length(data)
                if data[i][8]==data[j][8] && data[i][9]==data[j][9] && data[i][10]==data[j][10] && data[i][11]==data[j][11]
                    data[i][5]=data[i][5]+data[j][5]
                    data[i][6]=data[i][6]+data[j][6]
                    data[i][7]=data[i][7]+data[j][7]
                    push!(newdata,data[i])
                end
            end
        end
        data=newdata
    end

    #HyDe analysis
    data_length=length(data)
    for n in 1:data_length
        aabb_obs=data[n][5]
        abab_obs=data[n][6]
        abba_obs=data[n][7]
        mapname_t1=data[n][8]
        mapname_t2=data[n][9]
        mapname_t3=data[n][10]
        mapname_t4=data[n][11]

        p9 = (abba_obs+0.05)/num_obs #ABBA 
        p7 = (abab_obs+0.05)/num_obs #ABAB
        p4 = (aabb_obs+0.05)/num_obs #AABB

        #number of occurrences
        if havemap
            t1count=0
            t2count=0
            t3count=0
            t4count=0
            for i in 1:map_nrows
                if     mapd[i,2]==mapname_t1 t1count+=1
                elseif mapd[i,2]==mapname_t2 t2count+=1
                elseif mapd[i,2]==mapname_t3 t3count+=1                    
                elseif mapd[i,2]==mapname_t4 t4count+=1
                end
            end         
        else
            t1count=1
            t2count=1
            t3count=1
            t4count=1
        end                       
        
        avg_obs=num_obs/(t1count*t2count*t3count*t4count)

        #z-score
        obs_invp1 = avg_obs * (p9 - p7)
        obs_invp2 = avg_obs * (p4 - p7)
        obs_var_invp1 = avg_obs * p9 * (1 - p9) + avg_obs * p7 * (1 - p7) + 2 * avg_obs * p9 * p7
        obs_var_invp2 = avg_obs * p4 * (1 - p4) + avg_obs * p7 * (1 - p7) + 2 * avg_obs * p4 * p7
        obs_cov_invp1_invp2 = -1 * avg_obs * p9 * p4 + avg_obs * p9 * p7 + avg_obs * p7 * p4 + avg_obs * p7 * (1 - p7)
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

        #push to res
        if display_all
            if p_val<(alpha/3)
                ast="*"
                push!(res,[mapname_t1,mapname_t2,mapname_t3,mapname_t4,
                        aabb_obs,abab_obs,abba_obs,gamma,GH_ts,p_val,ast])
            else
                ast=""
                push!(res,[mapname_t1,mapname_t2,mapname_t3,mapname_t4,
                        aabb_obs,abab_obs,abba_obs,gamma,GH_ts,p_val,ast])
            end
        else
            if p_val<(alpha/3)
                ast="*"
                push!(res,[mapname_t1,mapname_t2,mapname_t3,mapname_t4,
                        aabb_obs,abab_obs,abba_obs,gamma,GH_ts,p_val,ast])
            end
        end
    end

    #prepare to print results in a cleaner way using DataFrame
    df=DataFrame(outgroup=String[],
                P1=String[],
                Hybrid=String[],
                P2=String[],
                AABB=Int[],
                ABAB=Int[],
                ABBA=Int[],
                Gamma=Float64[],
                Zscore=Float64[],
                Pvalue=Float64[],
                significance=String[])
    for result in res
        push!(df, result)
    end
    
    #write csv
    if isempty(filename) filename="HyDe-out" end
    if (writecsv)
        CSV.write("$filename.csv",df) 
        println("The results are stored as $filename.csv in the working directory.")    
    end

    println("Tip: if neccessary, use function showallDF(df) to see all the rows.")

    return df
    
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
