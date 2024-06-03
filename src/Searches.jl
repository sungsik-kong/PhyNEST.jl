const Move2Int = Dict{Symbol,Int}(:add=>1,:MVorigin=>2,:MVtarget=>3,:CHdir=>4,:delete=>5,:nni=>6)


# ProbST and add an optino to change starting points at each run
# add feature to ignore topologies when no root position presents
#seed (default 0 to get it from the clock): seed to replicate a given search
#if any error occurred, file .err provides information (seed) to reproduce the error.

"""
    preNNI(net::HybridNetwork,outgroup::AbstractString)

Function to make an NNI modification on the topology so the search can be initiated in different starting 
points. This could be very useful if we have low confidence in the starting topology.
In brief, using the input network, it unroots the topology, conduct NNI modification, removes branch lengths 
and set gamma to 0.5, and root the topology using the ougroup provided. We assume this should be (almost always) successful.

### Input
- **`net`**        A tree/network in HybridNetwork\\
- **`outgroup`**   Outgroup taxa in a string\\
"""
function preNNI(net::HybridNetwork,outgroup::AbstractString,nfindedge::Int64)
    suc=false
    nStartT=readTopologyUpdate(writeTopologyLevel1(net))#unroot the net
    while (!suc) suc=NNIRepeat!(nStartT,nfindedge) end
    for e in nStartT.edge e.length=-1.0 end
    @suppress begin nStartT=readTopology(writeTopologyLevel1(nStartT,outgroup)) end #Roots the network
    return nStartT
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
    
    ###########################
    #      Preliminaries      #
    ###########################
    #recycling the framework from SNaQ
    steps = 0
    failures = 0
    movescount = zeros(Int,18)
    movesfail = zeros(Int,6)
    Nmov = zeros(Int,6) 
    stillmoves = true
    outgroup_rooted=true #for topology checking during the search

    current_topology=starting_topology #setting N0
    res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera) #get optimization result (L0) and update bl and gamma on N0
    current_t_clikelihood=res.minimum #setting L0
    new_topology=deepcopy(current_topology) #copy N0
    
    ###########################
    #        Algorithm        #
    ###########################    
    while(steps<maximum_number_of_steps && failures<maximum_number_of_failures && stillmoves)
        steps+=1 #tracking number of steps in the space
        calculateNmov!(new_topology,Nmov) #function to calculate Nmov (i.e., number max of tries per move2int [0,0,0,0,0,0]=>[122,25,25,4,10000,42])
        move = whichMove(new_topology,hmax,movesfail,Nmov) #selects which move to take
        if move != :none
            #accepted=false
            new_topology=@suppress begin readTopologyUpdate(writeTopologyLevel1(new_topology,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) #unrooted, make modification on newT accroding to move, and get unrooted N1
            new_topology=@suppress begin readTopology(writeTopologyLevel1(new_topology,outgroup)) end #Roots the network using the outgroup, get N1
            #outgroup_rooted=correct_outgroup(new_topology,outgroup) #sometimes, the proposed topology cannot be rooted using the outgroup. Check if this happend.
            try rootatnode!(new_topology,outgroup)
            catch 
                outgroup_rooted=false
            end
            
            if(outgroup_rooted) #if the topology indeed rooted with the outgroup
                res,new_topology=do_optimization(new_topology,p,number_of_itera=number_of_itera) #get optimization result (L1) and update bl and gamma on N1
                new_t_clikelihood=res.minimum #setting L1
            else
                new_t_clikelihood=NaN #if the topology cannot be rooted using the outgroup, discard it.
            end
            
            #decide to accept N1 or not based on L1
            if new_t_clikelihood<current_t_clikelihood && !isnan(new_t_clikelihood)
                #accepted=true
                current_topology=new_topology
                current_t_clikelihood=new_t_clikelihood
                failures = 0
                movescount[Move2Int[move]+12] += 1
                movesfail = zeros(Int,6)#count of failed moves for current topology back to zero 
            else #accepted=false #not accepted
                new_topology=current_topology
                failures += 1
                movesfail[Move2Int[move]] += 1
            end
            #println(current_t_clikelihood)
        else
            #no more moves can proposed, so search ends
            stillmoves=false
        end
    end
    
    ###########################
    #         Results         #
    ###########################    
    str = "The search terminated at step $steps and at $(failures)th consecutive failures."
    str *="\nSummary of each move:"
    str *="\nInsertion of reticulation edge: $(movescount[1]) proposed, $(movescount[13]) accepted."
    str *="\nTail move of reticulation edge: $(movescount[2]) proposed, $(movescount[14]) accepted."
    str *="\nHead move of reticulation edge: $(movescount[3]) proposed, $(movescount[15]) accepted."
    str *="\nChange the direction of reticulation edge: $(movescount[4]) proposed, $(movescount[16]) accepted."
    str *="\nDeletion of reticulation edge: $(movescount[5]) proposed, $(movescount[17]) accepted."
    str *="\nNNI: $(movescount[6]) proposed, $(movescount[18]) accepted."
    str *="\nOn the current topology, $(sum(movesfail)) moves were made, including $(sum(movesfail)-failures) unsuccessful moves."
    if steps==maximum_number_of_steps
        str *="\nTerminated because it reached the maximum number of steps (current maximum_number_of_steps=$maximum_number_of_steps)."
        str *="\nIt is recommended to increase the number of maximum_number_of_steps and rerun the analysis.\n"
    elseif failures==maximum_number_of_failures
        str *="\nTerminated because it reached the maximum number of failures (current maximum_number_of_failures=$maximum_number_of_failures).\n"
    else 
        str *="\nTerminated although it neither reached the maximum number of steps or failures,"
        str *="\npossibly because there was no more move to make.\n"
    end
    print(str)
    if (write_log) writeflush(logfile,str) end

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

    ###########################
    #      Preliminaries      #
    ###########################
    steps=0
    movescount=zeros(Int,18)
    movesfail=zeros(Int,6)
    Nmov=zeros(Int,6)
    stillmoves=true
    
    current_topology=deepcopy(starting_topology)#N0 topology
    res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera)#get L0 and N0 phylogeny
    new_topology=deepcopy(current_topology)#N0
    current_t_clikelihood=res.minimum#L0
    
    #####################
    #      Burn-in      #
    #####################
    while(length(BurnIn)<number_of_burn_in)
        steps+=1
        calculateNmov!(new_topology,Nmov) #function to calculate Nmov (i.e., number max of tries per move2int [0,0,0,0,0,0]=>[122,25,25,4,10000,42])
        move = whichMove(new_topology,hmax,movesfail,Nmov) #selects which move to make

        if move != :none
            new_topology=@suppress begin readTopologyUpdate(writeTopologyLevel1(new_topology,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) 
            new_topology=@suppress begin readTopology(writeTopologyLevel1(new_topology,outgroup)) end #get the topology of N1
            new_t_clikelihood=0.0
            try
                res,new_topology=do_optimization(new_topology,p,number_of_itera=number_of_itera) #get L1 and updated N1
                new_t_clikelihood=res.minimum
            catch
                new_t_clikelihood=current_t_clikelihood
            end 
            
            
            delta=abs(current_t_clikelihood-new_t_clikelihood) #get differences

            current_t_clikelihood=new_t_clikelihood #there is no accept/reject criteria, always accept to make a change 
            if !isnan(delta) && !iszero(delta) && !isinf(delta) push!(BurnIn,delta) end #unless delta is weird... we use it
        else
            stillmoves=false
        end
    end

    return maximum(BurnIn) #find the largest delta (i.e., jump in likelihood)
end

"""
    numfails(p::Phylip, probN::Float64) 
"""
function numfails(p::Phylip, probN::Float64) 
    n=p.numtaxa
    numfail=(log(probN)/log(1-(1/(n-2)))/n)
    numfail=round(Int64,numfail)
    return numfail
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
    
    ###########################
    #      Preliminaries      #
    ###########################
    #recycling the framework from SNaQ plus some specifics for sim annealing
    ci=0.0
    steps=0
    ktrees=[]
    failures=0
    stillmoves=true
    outgroup_rooted=false #for topology checking during the search
    Nmov=zeros(Int,6) 
    movesfail=zeros(Int,6)#count of failed moves for current topology
    movescount=zeros(Int,18)#1:6 number of times moved proposed, 7:12 number of times success move (no intersecting cycles, etc.), 13:18 accepted by loglik
    
    current_topology=starting_topology #setting N0
    res,current_topology=do_optimization(current_topology,p,number_of_itera=number_of_itera) #get optimization result (L0) and update bl and gamma on N0
    current_t_clikelihood=res.minimum #setting L0
    #push N0 into the list of ktrees
    current_topology_newick=writeTopologyLevel1(current_topology) 
    push!(ktrees,(current_t_clikelihood,current_topology_newick))
    
    new_topology=deepcopy(current_topology) #copy N0
    
    ###########################
    #         Burn-in         #
    ###########################
    str="Running $number_of_burn_in burn-ins...\n"
    print(str)
    if (write_log) writeflush(logfile,str) end

    u=burn_in(current_topology,p,outgroup,hmax,number_of_burn_in,number_of_itera) 
    
    str="Burn in complete. At most -log likelihood value can change $u at a single step.\n"
    print(str)
    if (write_log) writeflush(logfile,str) end
    
    l=current_t_clikelihood
    b=cons/(((1-alph)*p.numtaxa)+((alph)*(log.(l)/p.seqleng)))
    
    ###########################
    #        Algorithm        #
    ###########################    
    while(steps < maximum_number_of_steps && failures < maximum_number_of_failures && stillmoves)
        steps+=1
        i=steps-1
        calculateNmov!(new_topology,Nmov)# function to calculate Nmov, number max of tries per move2int
        move = whichMove(new_topology,hmax,movesfail,Nmov)
        if move != :none
            outgroup_rooted=false
            accepted=false
            
            new_topology=@suppress begin readTopologyUpdate(writeTopologyLevel1(new_topology,di=true)) end #change to unrooted, plus update branch lengths to 1.0, we change to unrooted because proposedTop! requires it
            proposedTop!(move,new_topology,true,steps,10, movescount,movesfail,false) 
            new_topology=@suppress begin readTopology(writeTopologyLevel1(new_topology,outgroup)) end #get the topology of N1
            
            outgroup_rooted=correct_outgroup(new_topology,outgroup) #make sure N1 is rootable by the outgroup
            #try rootatnode!(new_topology,outgroup)
            #catch 
            #    outgroup_rooted=false
            #end

            if(outgroup_rooted) #if N1 is good to continue
                res,updated_new_topology=do_optimization(new_topology,p) #find L1 and update parameters on N1
                new_t_clikelihood=res.minimum
            else
                new_t_clikelihood=Inf #otherwise, discard it
            end
            
            #simulated annealing way of deciding to move
            if new_t_clikelihood<current_t_clikelihood && !isnan(new_t_clikelihood)
                accepted=true
            elseif new_t_clikelihood>current_t_clikelihood && !isnan(new_t_clikelihood) #deciding to aceept worse N1 at some probability
                ci=u/(1+(i*b))
                prob=exp((-1*(new_t_clikelihood-current_t_clikelihood))/ci)
                
                acceptance=[true,false]
                w=[prob,1-prob]
                accepted=sample(acceptance, Weights(w))
            else #discard otherwise
                accepted=false
            end
    
            if(accepted)
                #set N1 and L1 as new N0 and L0 and
                current_topology = new_topology
                current_t_clikelihood=new_t_clikelihood
                
                #reset steps and collect top accepted k Ns onto the list
                if !(new_t_clikelihood in ktrees[1]) #the accepted N1 is not in ktrees
                    failures = 0
                    movesfail = zeros(Int,6) 
                    movescount[Move2Int[move]+12] += 1
                    new_topology_newick=@suppress begin writeTopologyLevel1(updated_new_topology) end
                    klength=length(ktrees)
                    
                    if klength<k && new_t_clikelihood<ktrees[klength][1]
                        push!(ktrees,(new_t_clikelihood,new_topology_newick))
                    elseif klength==k && new_t_clikelihood<ktrees[klength][1]
                        ktrees[k]=(new_t_clikelihood,new_topology_newick)
                    else 
                        continue
                    end

                    ktrees=sort!(ktrees, by = x -> x[1])
                else #N1 is accepted but already in ktrees
                    #then, we consider it as a failed move to prevent endless searching
                    failures += 1
                    movesfail[Move2Int[move]] += 1
                end
               
            else
                #if not accepted, undo the move
                new_topology = current_topology
                failures += 1
                movesfail[Move2Int[move]] += 1
            end
        else
            stillmoves=false;
        end
    end

    ###########################
    #         Results         #
    ###########################    
    ktrees=unique(ktrees)
    rank=0
    str = "Showing $k best topologies found in this run:"
    str *="\nRank\tComposite Likelihood         Network"
    for i in ktrees
        rank+=1
    str *="\n$rank\t$(round(i[1],digits=5))         $(i[2])"
    end
    if(rank<k) str *="\nThe total length of best trees is shorter than k." end
    str *="\nSpeciation times for some newicks may not have updated if estimates are weird (e.g., NaN)."
    str *= "\nThe search terminated at step $steps and at $(failures)th consecutive failures and (Cooling schedule)ci=$ci."    
    str *= "\nSummary of each move:"
    str *= "\nInsertion of reticulation edge: $(movescount[1]) proposed, $(movescount[13]) accepted."
    str *= "\nTail move of reticulation edge: $(movescount[2]) proposed, $(movescount[14]) accepted."
    str *= "\nHead move of reticulation edge: $(movescount[3]) proposed, $(movescount[15]) accepted."
    str *= "\nChange the direction of reticulation edge: $(movescount[4]) proposed, $(movescount[16]) accepted."
    str *= "\nDeletion of reticulation edge: $(movescount[5]) proposed, $(movescount[17]) accepted."
    str *= "\nNearest-neighbor interchange (NNI): $(movescount[6]) proposed, $(movescount[18]) accepted."
    str *= "\nOn the current topology, $(sum(movesfail)) moves were made, including $(sum(movesfail)-failures) unsuccessful moves."

    
     #print reasons for termination
    if steps==maximum_number_of_steps
        str *= "\nTerminated because it reached the maximum number of steps (current maximum_number_of_steps=$maximum_number_of_steps)."
        str *= "\nIt is recommended to increase the number of maximum_number_of_steps and rerun the analysis.\n"
    elseif failures==maximum_number_of_failures
        str *= "\nTerminated because it reached the maximum number of failures (current maximum_number_of_failures=$maximum_number_of_failures).\n"
    else 
        str *= "\nTerminated although it neither reached the maximum number of steps or failures, possibly because there was no more move to make.\n"
    end
    
    print(str)
    if (write_log) writeflush(logfile,str) end

    #Return the first elements in ktrees array
    mpl=ktrees[1][1]
    bestnet= @suppress begin readTopology(ktrees[1][2]) end
    return mpl,bestnet
end

"""
    initiate_search(starting_topology::HybridNetwork,p::Phylip,outgroup::String,
        hmax::Integer,
        maximum_number_of_steps::Integer,
        maximum_number_of_failures::Integer,
        number_of_itera::Integer,
        number_of_runs::Integer,
        nniST::Bool,
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
                        nniprob::Float64,
                        do_hill_climbing::Bool,
                        number_of_burn_in::Int64,
                        k::Integer,
                        cons::Float64,
                        alph::Float64,
                        nfindedge::Int64,
                        write_log::Bool,
                        logfile::IO)

    search_results=[]

    if (do_hill_climbing)
        ###########################
        #      Hill Climbing      #
        ###########################
        search_results = Distributed.pmap(1:number_of_runs) do i            
            str="\n($i/$number_of_runs) Searching for the best network using the hill climbing algorithm; $(Dates.format(Dates.now(),"yyyy-mm-dd at HH:MM:SS"))\n"
            print(str)
            if (write_log) writeflush(logfile,str) end

            #NNI strating topology
            probnni=rand(Uniform(0,1))
            if(probnni<nniprob)
                starting_topology=preNNI(starting_topology,outgroup,nfindedge) 
                str="Starting topology changed to $(writeTopologyLevel1(starting_topology)) using nearest-neighbor interchange (NNI)\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
            else
                str="Starting topology remains unchanged\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
            end

            #executing the algorithm
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
                if (write_log) writeflush(logfile,str) end
                
                push!(search_results,(current_t_clikelihood,current_topology))
                return current_t_clikelihood,current_topology
            catch err 
                str="($i/$number_of_runs) Terminated due to an error\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
                #throw(err)
                @error "ERROR: " exception=(err, catch_backtrace())
            end
        end
    else
        ###########################
        #      Sim Annealing      #
        ###########################
        search_results = Distributed.pmap(1:number_of_runs) do i
            #i+=1      
            str=("\n($i/$number_of_runs) Searching for the best network using the simulated annealing algorithm; $(Dates.format(Dates.now(),"yyyy-mm-dd at HH:MM:SS"))\n")
            print(str)
            if (write_log)
                write(logfile,str)
                flush(logfile)
            end

            #NNI strating topology            
            probnni=rand(Uniform(0,1))
            if(probnni<nniprob)
                starting_topology=preNNI(starting_topology,outgroup,nfindedge) 
                str="Starting topology changed to $(writeTopologyLevel1(starting_topology)) using nearest-neighbor interchange (NNI)\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
            else
                str="Starting topology remains unchanged\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
            end
            
            #executing the algorithm
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
                if (write_log) writeflush(logfile,str) end
                return current_t_clikelihood,current_topology
            catch err
                str="($i/$number_of_runs) Terminated due to an error\n"
                print(str)
                if (write_log) writeflush(logfile,str) end
                @error "ERROR: " exception=(err, catch_backtrace())
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
            nniST=false::Bool,
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
- `nniST (default=false)`                       Conducts NNI operation on the specified starting topology if set as `true`

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
                maximum_number_of_failures=150::Integer,
                probN=9.5e-45::Float64,
                number_of_itera=1000::Integer,
                number_of_runs=10::Integer,
                nniprob=0.75::Float64,
                do_hill_climbing=true::Bool,
                number_of_burn_in=25::Integer,
                k=10::Integer,
                cons=0.5::Float64,
                alph=0.5::Float64,
                nfindedge=10::Int64,
                filename=""::AbstractString
                )

    ###########################
    #      Preliminaries      #
    ###########################
    #current time
    now=Dates.now()

    #output filename
    if isempty(filename) filename="PhyNEST" end
    if(do_hill_climbing) filename=filename*("_hc") 
    else filename=filename*("_sa") end

    #open log file
    log=string(filename,".log")
    logfile=open(log,"w")
    
    #number of processors when using parallel computing
    #to decide whether we can write stuff in the logfile
    write_log=true
    num_processors=Distributed.nprocs()
    if num_processors!==1 write_log=false end

    #maximum number of failures for termination plus
    #which searching strategy algorithm selected
    if (do_hill_climbing) 
        maximum_number_of_failures=maximum_number_of_failures
        maxfailcomment="using maximum_number_of_failures provided"
        algorithm="Hill climbing"
    else 
        maximum_number_of_failures=numfails(p, probN) 
        maxfailcomment="computed using probN=$probN provided."
        algorithm="Simulated annealing"
    end

    ############################
    #      Printing Basics     #
    ############################
    str="PhyNEST: Phylogenetic Network Estimation using SiTe patterns"
    str*="\nAnalysis start: $(Dates.format(now, "yyyy-mm-dd at HH:MM:SS"))"
    str*="\nInput file: $(p.filename)"
    str*="\nNumber of sequences: $(p.numtaxa)"
    str*="\nSequence length: $(p.seqleng)"
    str*="\nStarting Topology: $(writeTopologyLevel1(starting_topology))"
    str*="\nOutgroup specified for rooting: $outgroup"
    str*="\nNumber of maximum reticulation(s): $hmax"
    str*="\nThe maximum number of iterations for each optimization: $number_of_itera"
    str*="\nSearch algorithm selected: $algorithm "
    if !(do_hill_climbing)
        str*="with alpha=$alph; c=$cons"
    end
    str*="\nThe maximum number of steps during search: $maximum_number_of_steps"
    str*="\nThe maximum number of consecutive failures: $maximum_number_of_failures ($maxfailcomment)"
    str*="\nOutput file store directory: $(pwd())"
    str*="\nThe number of processors for this analysis: $(num_processors)\n"
    printwriteflush(logfile,str)


    ########################
    #      Searching       #
    ########################
    search_results=initiate_search(starting_topology,p,outgroup,
                                    hmax,
                                    maximum_number_of_steps,
                                    maximum_number_of_failures,
                                    number_of_itera,
                                    number_of_runs,
                                    nniprob,
                                    do_hill_climbing,
                                    number_of_burn_in,
                                    k,
                                    cons,
                                    alph,
                                    nfindedge,
                                    write_log,
                                    logfile
                                    )

#    display(search_results)                        

    #########################
    #     Summarization     #
    #########################
    run=0    
    str="\nSummary:\n"
    str*=
    "Run\tComposite likelihood\tNetwork estimated\n"
        for i in search_results
            run+=1
            if i!==nothing
                clik=round(i[1],digits=1)    
                str *=
                "$run\t$clik\t\t$(writeTopologyLevel1(i[2]))\n"
            else
                str *=
                "$run\tRun $run was aborted due to an error\n"
            end
        end
        printwriteflush(logfile,str)
    
    deleteat!(search_results, findall(x->x==nothing, search_results))

    if isempty(search_results)
        return "None of the heuristic searches were successful. Check the search parameter settings and try again."
    else

    sorted_search_results=sort(search_results, by = first)
    best_topology=sorted_search_results[1][2]
    best_top=writeTopologyLevel1(best_topology)
    best_top_di=writeTopologyLevel1(best_topology,di=true)

    str="\nBest topology: $best_top)\n"
    str*="End;\n"
    printwriteflush(logfile,str)

    #############################
    #     Creating .out file    #
    #############################

    out=string(filename,".out")
    outfile=open(out,"w")
    write(outfile,"$best_top\n\n")
    write(outfile,best_top_di)
    flush(outfile)

    return best_topology
    end
end

function printwriteflush(logfile::IO,str::String)
    print(str)
    write(logfile,str)
    flush(logfile)
end

function writeflush(logfile::IO,str::String)
    write(logfile,str)
    flush(logfile)
end