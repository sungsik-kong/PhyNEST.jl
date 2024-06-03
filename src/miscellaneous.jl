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
    now=Dates.now()
    println("""
    Thank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.
    Please report bugs or make suggestions to https://groups.google.com/g/phynest-users.
    If you conduct an analysis using PhyNEST, please cite:
    Sungsik Kong, David Swofford, and Laura Kubatko (2022) Inference of Phylogenetic Networks from Sequence Data using Composite Likelihood.
    Preprint available online on BioRxiv at https://doi.org/10.1101/2022.11.14.516468.
    time stamp: $(Dates.format(now, "yyyy-mm-dd at HH:MM:SS"))
    """)
end

"""
    get_average(i::Array)

Compute average of the values in an array. 
"""
function get_average(i::Array) return sum(i)/length(i) end












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
- `writecsv (default=false)` If `true`, the result is exported and stored in `.csv` file in the working directory
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
function showallDF(df::DataFrame) CSV.show(df,allrows=true) end

"""
    HyDe(p::Phylip, outgroup::AbstractString; 
        alpha=0.05::Float64, 
        display_all=true::Bool, #filter
        map=""::AbstractString,
        writecsv=false::Bool, 
        filename=""::AbstractString


Conducts HyDe: Hybrid Detection using Phylogenetic invariants. See Blischak et al., (2018) (https://doi.org/10.1093/sysbio/syy023) and the manual (https://hybridization-detection.readthedocs.io) for more information. 
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
2×11 DataFrame
 Row │ outgroup  P1      Hybrid  P2      AABB   ABAB   ABBA   Gamma     Zscore   Pvalue   significance
     │ String    String  String  String  Int64  Int64  Int64  Float64   Float64  Float64  String
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 5         3       2       1        8005   1991   8057  0.502152  47.6571      0.0  *
   2 │ 5         1       2       3        8057   1991   8005  0.497848  47.6571      0.0  *

julia> HyDe(p,"5",map="map.txt",display_all=false)
Map file [map.txt] provided.
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

    display(mapdict)

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
                                t4]
                        )
                end
            end
        end
    end

    #summarize and add up according to the map file
    if havemap && length(unique(keys(mapdict)))<p.numtaxa
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

    return df
    
end

"""
    function cct(pvals::AbstractArray; weights=Float64[]::Array)

A function to perform the Cauchy combination test. It takes a list of
p-values and a list of weights and return the global p-value. Example below shows how the p-values
computed from the function `HyDe` can be directly used for `cct`
        
See [https://github.com/rhaque62/pyghdet] for more information. This function is the
`PhyNEST` implementation of the Cauchy combination test proposed in Haque and Kubatko (2023)
[https://www.biorxiv.org/content/biorxiv/early/2023/02/27/2023.02.24.529943.full.pdf],
a method that consider hybridization and coalescence in a unified framework that can
detect whether there are any hybrid species in a given set of species. Based on this global
test of hybridization, one can decide whether a tree or network analysis is appropriate for 
a given data set.
"""
function cct(pvals::Array; weights=Float64[]::Array,
            lower_bound=1e-17::Float64, upper_bound=0.99::Float64)

    #check if there is any non-numeric values in the pvals
    all_float=Float64 in typeof.(pvals)
    all_float || error("The individual tests produced p-values containing non-numeric character! Failed to test the global null hypothesis")
    
    #check if all element in pv_arr < 1 and > 0
    p_btw_zero_and_one=all(>=(0.0), pvals) && all(<=(1.0), pvals)
    p_btw_zero_and_one || error("All the individual p-values must be between 0 and 1! Failed to test the global null hypothesis")
    
    # making the p-value ready for cct
    pv=Float64[]
    for p in pvals
        if p==0 push!(pv,lower_bound)
        elseif p==1 push!(pv,upper_bound)
        else push!(pv,p)
        end
    end
    
    #check the weights length=pval length; all weight > 0
    if isempty(weights)
        weights=zeros(length(pvals))
        w=1/length(pvals)
        fill!(weights, w)
    end
    #println(weights)
    all(>=(0), weights) && all(<=(1), weights) ||  #&& sum(weights)==1 
        error("There are some problems with the weights. Check again. (weights=$weights)")

    #cct computation
    small_p=Float64[]
    small_w=Float64[]
    large_p=Float64[]
    large_w=Float64[]

    i=1
    for item in pv
        if item < (lower_bound*10)
            push!(small_p,item)
            push!(small_w,weights[i])
            i+=1
        else
            push!(large_p,item)
            push!(large_w,weights[i])
        end
    end

    if length(small_p)==0
        cct_stat=[]
        j=1
        for item in large_p
            np=large_w[j]*tan((0.5-item)*pi)
            push!(cct_stat,np)
            j+=1
        end
        cct_stat=sum(cct_stat)
    else
        cct_small=Float64[]
        cct_large=Float64[]
        k=1
        for item in small_p
            np=(weights[k]/item)/pi
            push!(cct_small,np)
            k=+1
        end
        cct_small=sum(cct_small)

        l=1
        for item in large_p
            np=large_w[l]*tan(0.5-item)*pi
            push!(cct_large,np)
            l+=1
        end

        cct_large=sum(cct_large)
        cct_stat=cct_small+cct_large
    end

    #calculate the pv for the global test
    D=Distributions.Cauchy()
    if cct_stat>1e15
        pval=(1/cct_stat)/pi
    else
        pval=1-Distributions.cdf(D,cct_stat)
    end

    return pval
end

"""
    function cmc(pvals::Array)

A function to perform the CMC test. It takes a list of p-values and return the global p-value.
"""
function cmc(pvals::Array)
    p_min = min(1, length(pvals)*findmin(pvals)[1])
    p_cmc=cct([cct(pvals),p_min])
    
    return p_cmc
end



"""
    function LRT
"""

function LRT(p::Phylip)
    s2=readTopology("(O,((P1,H),P2));") #gamma=0
    s=readTopology("(O,((P2,(H)#H1:::0.51),(P1,#H1:::0.49)));") #gamma!==0
    #display(p)
    #display(s2)
    #display(s)

    sups2=do_optimization(s2,p)
    sups=do_optimization(s,p)
    
    sups2=(sups2[1].minimum)
    sups=(sups[1].minimum)

    lambda=-2*(sups2-sups)

    println(lambda)
end
