const Move2Int = Dict{Symbol,Int}(:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5,:nni=>6)

###---heuristics using hill climbing algorithm---###
"""
    hill_climbing(starting_topology::HybridNetwork,p::Phylip,outgroup::String;
        hmax::Integer,
        maximum_number_of_steps::Integer,
        maximum_number_of_failures::Integer,
        number_of_itera::Integer)

Executes hill climbing algorithm for the network space search.
"""
function hill_climbing(starting_topology::HybridNetwork,p::Phylip,outgroup::String,
                        hmax::Integer,
                        maximum_number_of_steps::Integer,
                        maximum_number_of_failures::Integer,
                        number_of_itera::Integer,
                        write_log::Bool,
                        logfile::IO)
#set some necessary stuff for search
steps = 0
failures = 0
movescount = zeros(Int,18)#1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
movesfail = zeros(Int,6)#count of failed moves for current topology
Nmov = zeros(Int,6) 
stillmoves = true

#set current topology, copy it and set as new topology, and clikelihood for the current topology
#current_topology=readTopology(writeTopologyLevel1(starting_topology))
current_topology=starting_topology

#if probST<1.0 && rand() < 1-probST # modify starting tree by a nni move
#    suc = NNIRepeat!(current_topology,10); #will try 10 attempts to do an nni move, if set to 1, hard to find it depending on currT
#end


res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera,update_parameters=true)
current_t_clikelihood=res.minimum
new_topology=deepcopy(current_topology)

#begin search
while(steps<maximum_number_of_steps && failures<maximum_number_of_failures && stillmoves)
    steps+=1 #tracking number of topologies evaluated 
    calculateNmov!(new_topology,Nmov)# function to calculate Nmov, number max of tries per move2int [0,0,0,0,0,0]=>[122,25,25,4,10000,42]
    move = whichMove(new_topology,hmax,movesfail,Nmov) #selects which move to take
    if move != :none
        outgroup_rooted=false
        accepted=false
        
        new_topology=@suppress begin readTopologyUpdate(writeTopologyLevel1(new_topology,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
        proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) #unrooted, make modification on newT accroding to move
        new_topology=@suppress begin readTopology(writeTopologyLevel1(new_topology,outgroup)) end#Roots the network
        outgroup_rooted=correct_outgroup(new_topology,outgroup)
             
        if(outgroup_rooted)
            res,new_topology=do_optimization(new_topology,p,number_of_itera=number_of_itera,update_parameters=true)
            new_t_clikelihood=res.minimum
        else
            new_t_clikelihood=Inf
        end
        
        if(new_t_clikelihood<current_t_clikelihood)
            accepted=true
        else
            accepted=false
        end

        if(accepted)
            current_topology=new_topology
            current_t_clikelihood=new_t_clikelihood
            failures = 0
            movescount[Move2Int[move]+12] += 1
            movesfail = zeros(Int,6)#count of failed moves for current topology back to zero 
        else
            new_topology=current_topology
            failures += 1
            movesfail[Move2Int[move]] += 1
        end
    else
        stillmoves=false
    end
end

#determine the reason for termination
termination=0
if steps==maximum_number_of_steps termination=1
elseif failures==maximum_number_of_failures termination=2
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
"Terminated because it reached the maximum number of steps (current maximum_number_of_steps=$maximum_number_of_steps). 
It is recommended to increase the number of maximum_number_of_steps and rerun the analysis.\n"
elseif termination==2 
    str *= 
"Terminated because it reached the maximum number of failures (current maximum_number_of_failures=$maximum_number_of_failures).\n"
else 
    str *= 
"Terminated although it neither reached the maximum number of steps or failures,
possibly because there was no more move to make.\n"
end

print(str)
if(write_log)
#    print(str)
    write(logfile,str)
    flush(logfile)
end

return current_t_clikelihood,current_topology

end

"""
    burn_in(starting_topology::HybridNetwork, p::Phylip, outgroup::String, 
        hmax::Integer, number_of_burn_in::Int64, number_of_itera::Integer)

Burn in procedure that needs to be conducted in prior to the simulated annealing network search.
"""
function burn_in(starting_topology::HybridNetwork, p::Phylip, outgroup::String, 
                    hmax::Integer, number_of_burn_in::Int64, number_of_itera::Integer)
    BurnIn=Float64[]

    steps = 0
    movescount = zeros(Int,18)
    movesfail = zeros(Int,6)
    Nmov = zeros(Int,6)
    stillmoves = true
    

    current_topology=deepcopy(starting_topology)
    res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera,update_parameters=false)
    current_t_clikelihood=res.minimum
    new_topology=deepcopy(current_topology)

    while(length(BurnIn)<number_of_burn_in)
        steps+=1

        calculateNmov!(new_topology,Nmov)# function to calculate Nmov, number max of tries per move2int
        move = whichMove(new_topology,hmax,movesfail,Nmov)

        if move != :none
            new_topology=readTopologyUpdate(writeTopologyLevel1(new_topology)) #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) #unrooted, make modification on newT accroding to move
            @suppress begin new_topology=readTopology(writeTopologyLevel1(new_topology,outgroup)) end #Roots the network
            res,new_topology=do_optimization(new_topology,p,number_of_itera=number_of_itera,update_parameters=false)
            new_t_clikelihood=res.minimum
            
            delta=abs(current_t_clikelihood-new_t_clikelihood)

            current_t_clikelihood=new_t_clikelihood
            if !isnan(delta) && !iszero(delta) && !isinf(delta)
                push!(BurnIn,delta)
            end
        else
            stillmoves=false
        end
    end
    
    return findmax(BurnIn)[1] 
end

"""
    simulated_annealing(starting_topology::HybridNetwork, p::Phylip, outgroup::String, 
        hmax::Integer, 
        number_of_burn_in::Int64,
        maximum_number_of_steps::Integer,
        maximum_number_of_failures::Integer,
        number_of_itera::Integer,
        k::Integer,
        cons::Float64,
        alph::Float64
        )

Executes simulated algorithm for the network space search.
"""
function simulated_annealing(starting_topology::HybridNetwork, p::Phylip, outgroup::String, 
                                hmax::Integer, 
                                number_of_burn_in::Int64,
                                maximum_number_of_steps::Integer,
                                maximum_number_of_failures::Integer,
                                number_of_itera::Integer,
                                k::Integer,
                                cons::Float64,
                                alph::Float64,
                                write_log::Bool,
                                logfile::IO)
    
    #set some necessary stuff for search
    steps = 0
    failures = 0
    movescount = zeros(Int,18)#1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
    movesfail = zeros(Int,6)#count of failed moves for current topology
    Nmov = zeros(Int,6) 
    stillmoves = true
    #--new stuffs that are not present in hill climbing
    ktrees=[]
    ci=0.0

    #dealing with the starting tree
    current_topology=starting_topology
    #current_topology=readTopology(writeTopologyLevel1(starting_topology))
    res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera,update_parameters=true)
    current_t_clikelihood=res.minimum
    current_topology_newick=writeTopologyLevel1(current_topology)
    push!(ktrees,(current_t_clikelihood,current_topology_newick))
    #println(ktrees)

    new_topology=deepcopy(current_topology)

    #Burn in
    str="Running $number_of_burn_in burn ins...\n"
    print(str)
    if (write_log)
        write(logfile,str)
        flush(logfile)
    end
    u=burn_in(current_topology,p,outgroup,hmax,number_of_burn_in,number_of_itera) 
    str="Burn in complete. At most -log likelihood value can change $u at a single step.\n"
    print(str)
    if (write_log)
        write(logfile,str)
        flush(logfile)
    end
    l=current_t_clikelihood
    b=cons/(((1-alph)*p.numtaxa)+((alph)*(log.(l)/p.seqleng)))
    

    while(steps < maximum_number_of_steps && failures < maximum_number_of_failures && stillmoves)
        steps+=1
        i=steps-1
        calculateNmov!(new_topology,Nmov)# function to calculate Nmov, number max of tries per move2int
        move = whichMove(new_topology,hmax,movesfail,Nmov)
        if move != :none
            outgroup_rooted=false
            accepted=false
            #println(new_topology)
            #new_topology=readTopologyUpdate(writeTopologyLevel1(new_topology))
            new_topology=@suppress begin readTopologyUpdate(writeTopologyLevel1(new_topology,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) 
            new_topology=@suppress begin readTopology(writeTopologyLevel1(new_topology,outgroup)) end
            outgroup_rooted=correct_outgroup(new_topology,outgroup)
            #println("hi")
            if(outgroup_rooted)
                res,updated_new_topology=do_optimization(new_topology,p,update_parameters=true)
                new_t_clikelihood=res.minimum
            else
            new_t_clikelihood=Inf
            end
            
            if new_t_clikelihood<=current_t_clikelihood && !isnan(new_t_clikelihood) && !isnan(current_t_clikelihood)
                accepted=true
            elseif new_t_clikelihood>current_t_clikelihood && !isnan(new_t_clikelihood) && !isnan(current_t_clikelihood)
                ci=u/(1+(i*b))
                prob=exp((-1*(new_t_clikelihood-current_t_clikelihood))/ci)
                items=[true,false]
                weigh=[prob,1-prob]
                accepted=sample(items, Weights(weigh))
            else
                accepted=false
            end

            if(accepted)
                new_topology_newick=writeTopologyLevel1(updated_new_topology)
                ktrees=sort!(ktrees, by = x -> x[1])
                #println(ktrees)
                flag=false
                for e in ktrees
                    if (new_t_clikelihood in e[1])==false
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
                    #if steps==1
                    #    push!(ktrees,(current_t_clikelihood,current_topology_newick))
                        #println(ktrees)
                    #else
                        if flag==false
                            push!(ktrees,(new_t_clikelihood,new_topology_newick))
                            #println(ktrees)
                        end
                    #end
                elseif length(ktrees)==k && flag==false
                    if new_t_clikelihood<ktrees[k][1]
                        ktrees[k]=(new_t_clikelihood,new_topology_newick)
                        #println(ktrees)
                    end
                end
                current_topology = new_topology
                current_t_clikelihood=new_t_clikelihood
            else
                new_topology = current_topology
                failures += 1
                movesfail[Move2Int[move]] += 1
            end
        else
            stillmoves=false;
        end
    end

    ktrees=sort!(ktrees, by = x -> x[1])

    rank=0
    str =
"Showing $k best topologies found in this run:
Rank   Composite Likelihood    Network\n"
    for i in ktrees
        rank+=1
    str *=
"$rank\t$(round(i[1],digits=5))         $(i[2])\n"
    end
    
    if(rank<k) str *="The total length of best trees can be shorter than k.\n" end
    str *="Speciation times for some newicks may not have updated if estimates are weird (e.g., NaN).\n"
    print(str)
    if(write_log)
        write(logfile,str)
        flush(logfile)
    end
    mpl=ktrees[1][1]
    bestnet=readTopology(ktrees[1][2])
    
    #determine the reason for termination
     termination=0
     if steps==maximum_number_of_steps termination=1
     elseif failures==maximum_number_of_failures termination=2
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
"Terminated because it reached the maximum number of steps (current maximum_number_of_steps=$maximum_number_of_steps). 
It is recommended to increase the number of maximum_number_of_steps and rerun the analysis.\n"
    elseif termination==2 
        str *= 
"Terminated because it reached the maximum number of failures (current maximum_number_of_failures=$maximum_number_of_failures).\n"
    else 
        str *= 
"Terminated although it neither reached the maximum number of steps or failures,
possibly because there was no more move to make.\n"
    end
    
    print(str)
    if(write_log)    
        write(logfile,str)
        flush(logfile)
    end
    
    return mpl,bestnet

end

"""
    initiate_search(starting_topology::HybridNetwork,p::Phylip,outgroup::String,
        hmax::Integer,
        maximum_number_of_steps::Integer,
        maximum_number_of_failures::Integer,
        number_of_itera::Integer,
        number_of_runs::Integer,
        do_hill_climbing::Bool,
        number_of_burn_in::Int64,
        k::Integer,
        cons::Float64,
        alph::Float64)

Depending on the setting provided in phyne!, either conducts hill climbing or simulated annealing searching.
"""
function initiate_search(starting_topology::HybridNetwork,p::Phylip,outgroup::String,
    hmax::Integer,
    maximum_number_of_steps::Integer,
    maximum_number_of_failures::Integer,
    number_of_itera::Integer,
    number_of_runs::Integer,
    do_hill_climbing::Bool,
    number_of_burn_in::Int64,
    k::Integer,
    cons::Float64,
    alph::Float64,
    write_log::Bool,
    logfile::IO)

    search_results=[]
    #i=1

    if (do_hill_climbing)
        #while i < number_of_runs
        #for i in 1:number_of_runs
        search_results = Distributed.pmap(1:number_of_runs) do i
            str="\n($i/$number_of_runs) Searching for the best network using the hill climbing algorithm, $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n"
            print(str)
            if (write_log)
                write(logfile,str)
                flush(logfile)
            end
            try
                current_t_clikelihood,current_topology=hill_climbing(starting_topology,p,outgroup,
                                                                    hmax,
                                                                    maximum_number_of_steps,
                                                                    maximum_number_of_failures,
                                                                    number_of_itera,
                                                                    write_log,
                                                                    logfile)
                
                str="($i/$number_of_runs) Estimated topology in this run: $(writeTopologyLevel1(current_topology))\n"
                str*="($i/$number_of_runs) Composite likelihood of the estimated topology in this run: $(current_t_clikelihood)\n"
                print(str)
                if (write_log)
                    write(logfile,str)
                    flush(logfile)
                end
                push!(search_results,(current_t_clikelihood,current_topology))
                return current_t_clikelihood,current_topology
            catch 
                str="($i/$number_of_runs) Terminated due to error\n"
                print(str)
                if (write_log)
                    write(logfile,str)
                    flush(logfile)
                end
            end
        end
        
        #println("best: $best")
        #i+=1
    else
        #for i in 1:number_of_runs
        search_results = Distributed.pmap(1:number_of_runs) do i
            #i+=1      
            str=("\n($i/$number_of_runs) Searching for the best network using the simulated annealing algorithm, $(Dates.format(Dates.now(), "yyyy-mm-dd H:M:S.s"))\n")
            print(str)
            if (write_log)
                write(logfile,str)
                flush(logfile)
            end
            try 
                current_t_clikelihood,current_topology=simulated_annealing(starting_topology, p, outgroup, 
                            hmax,
                            number_of_burn_in,
                            maximum_number_of_steps,
                            maximum_number_of_failures,
                            number_of_itera,
                            k,
                            cons,
                            alph,
                            write_log,
                            logfile)
                push!(search_results,(current_t_clikelihood,current_topology))               
                str=("($i/$number_of_runs) Estimated topology in this run: $(writeTopologyLevel1(current_topology))\n")
                str*=("($i/$number_of_runs) Composite likelihood of the estimated topology in this run: $(current_t_clikelihood)\n")
                print(str)
                if (write_log)
                    write(logfile,str)
                    flush(logfile)
                end
                return current_t_clikelihood,current_topology
            catch 
                str="($i/$number_of_runs) Terminated due to error\n"
                print(str)
                if (write_log)
                    write(logfile,str)
                    flush(logfile)
                end
            end
        end
    end

    return search_results
end




















































































































































"""
    phyne!(starting_topology::HybridNetwork,p::Phylip,outgroup::String;
            hmax=1::Integer,
            maximum_number_of_steps=250000::Integer,
            maximum_number_of_failures=100::Integer,
            number_of_itera=1000::Integer,
            number_of_runs=10::Integer,
            do_hill_climbing=true::Bool,
            number_of_burn_in=25::Integer,
            k=10::Integer,
            cons=0.5::Float64,
            alph=0.5::Float64,
            filename=""::AbstractString)

`phyne!` is fine. `phyne!` executes function `initiate_search(args)`.

Estimate the species network (or tree if `hmax=0`) using maximum composite likelihood. The search begins from the `starting_topology`
which can be either estimated or randomly generated. Starting topology can be either tree or a network with `<=hamx`. Outgroup taxon 
must be specified to root the network. 

## Mandatory arguments
- `starting_topology`       Starting topology in `HybridNetwork` object created using the function `readTopology` 
- `p`                       Sequence alignment parsed in `Phylip` object using the function `readPhylip`
- `outgroup`                Name of the outgroup taxa

## Optional arguments
### Generic
- `hmax (dafault=1)`                            Maximum number of hybridizations to be included in the final network
- `do_hill_climbing (default=true)`             When `=true`, network is searched using hill climbing and when `=false`, it searches using simulated annealing
- `number_of_runs (default=10)`                 Number of independent runs

### Network space search
- `maximum_number_of_steps (default=250000)`    Maximum number of steps before the search terminates
- `maximum_number_of_failures (default=100)`    Maximum number of consecutive failures (i.e., rejecting the proposed topology) before the search teminates

### Optimization 
- `number_of_itera (default=1000)`

### For simulated annealing
- `number_of_burn_in (default=25)`              
- `k (default=10)`                              Specifies the number of best networks to be stored at each run
- `cons (default=0.5)`                          
- `alph (default=0.5)`

### Output
- `filename (default="")` Specifies the name of the output file. If unspecified, it will use `PhyNEST_hc` or `PhyNEST_sa` depending on the heuristic method applied.
- The best network estimated throughout the entire runs is written in the extended Newick format and stored in file with an extension `.out`. 
First line in the file is readable by `PhyloNetworks` that can be visualized using `PhyloPlots`. Third (last) line is the identical network
but readable by `DendroScope`. 
- Specifics of the network search from each independent run is stored in the file with an extension `.log`. 
It summarizes how many of each network 'moves' were made during the search, the reason for terminating the search 
(e.g., researched maximum number of steps, reached maximum number of consecutive failures, or else), and 
the final network estimated and its composite likelihood. When `do_hill_climbing=true`, a single best network is selected for each run,
but when `do_hill_climbing=false` it prints `k` best networks visited.
"""
function phyne!(starting_topology::HybridNetwork,p::Phylip,outgroup::String;
    hmax=1::Integer,
    maximum_number_of_steps=250000::Integer,
    maximum_number_of_failures=100::Integer,
    number_of_itera=1000::Integer,
    number_of_runs=10::Integer,
    do_hill_climbing=true::Bool,
    number_of_burn_in=25::Integer,
    k=10::Integer,
    cons=0.5::Float64,
    alph=0.5::Float64,
    filename=""::AbstractString
    )

    #current time
    t=Dates.now()

    #filename
    if isempty(filename) filename="PhyNEST" end
    if(do_hill_climbing) filename=filename*("_hc") 
    else filename=filename*("_sa") end

    #open log file
    log=string(filename,".log")
    logfile=open(log,"w")
    typeof(logfile)
    #if(display) print(str) end
    
    #number of processors when using parallel computing
    #to decide whether we can write stuff in the logfile
    num_processors=Distributed.nprocs()
    write_log=true
    if num_processors!==1 write_log=false end
    
    #which searching strategy algorithm selected
    if (do_hill_climbing) algorithm="Hill climbing"
    else algorithm="Simulated annealing" end

str="PhyNEST: Phylogenetic Network Estimation using SiTe patterns
Analysis start: $(Dates.format(t, "yyyy-mm-dd at HH:MM:SS"))
Input sequence file: $(p.filename)
Number of sequences: $(p.numtaxa) 
Sequence length: $(p.seqleng)
Starting Topology: $(writeTopologyLevel1(starting_topology))
Outgroup specified for rooting: $outgroup
Number of maximum reticulation(s): $hmax
The maximum number of iterations for each optimization: $number_of_itera
Search algorithm selected: $algorithm
The maximum number of steps during search: $maximum_number_of_steps
Output file store path: $(pwd())
The number of processors for this analysis: $(num_processors)\n"
    print(str)
    write(logfile,str)
    flush(logfile)
    @debug "[$(Dates.now())] Network searching using PhyNEST is starting"
    
    search_results=initiate_search(starting_topology,p,outgroup,
                        hmax,
                        maximum_number_of_steps,
                        maximum_number_of_failures,
                        number_of_itera,
                        number_of_runs,
                        do_hill_climbing,
                        number_of_burn_in,
                        k,
                        cons,
                        alph,
                        write_log,
                        logfile
                        )

    run=0    
    str="\n-----Summary of the networks found\n"
    str*=
    "Run   Composite Likelihood    Network\n"
        for i in search_results
            run+=1
        str *=
    "$run\t$(round(i[1],digits=5))         $(writeTopologyLevel1(i[2]))\n"
        end
    print(str)
    write(logfile,str)
    flush(logfile)

    sorted_search_results=sort(search_results, by = first)
    @debug "[$(Dates.now())] The search results from $(number_of_runs) runs are sorted by their composite likelihood"
    best_topology=sorted_search_results[1][2]
    @debug "[$(Dates.now())] The 'best' topology is selected as a final product stored"

str="\n-----end of analysis\n"
    print(str)

    best_top=writeTopologyLevel1(best_topology)
    best_top_di=writeTopologyLevel1(best_topology,di=true)

str*="Best topology: $best_top)"
    write(logfile,str)
    flush(logfile)

    out=string(filename,".out")
    outfile=open(out,"w")
    write(outfile,"$best_top\n\n")
    write(outfile,best_top_di)
    flush(outfile)

    return best_topology
end


# ProbST and add an optino to change starting points at each run
# add feature to ignore topologies when no root position presents
#seed (default 0 to get it from the clock): seed to replicate a given search
#if any error occurred, file .err provides information (seed) to reproduce the error.
