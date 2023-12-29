#written by Sungsik Kong 2021-2022
"""
    do_optimization(net::HybridNetwork, p::Phylip;
                    lower_bound=0.00001::Float64,
                    factor=2.0::Float64,
                    tolerance=0.01::Float64,
                    number_of_itera=1000::Int64)

Optimizes the composite likelihood of the network topology given the data.


"""
function do_optimization(topology::HybridNetwork, p::Phylip;
                        lower_bound=0.00001::Float64,
                        factor=2.0::Float64,
                        tolerance=0.01::Float64,
                        number_of_itera=1000::Int64)
    
    net=deepcopy(topology)

    #---check topology---#
    #weather or not gamma values are specified, if not, arbitrarily insert 1/2*h.
    for e in net.edge
        if e.hybrid && !(0<e.gamma<1)
            e.gamma=1/(2*net.numHybrids)
        end
    end
    
    #---preamble---#
    parameter_sets=Vector{Array{Float64}}()#Array[]   
    #decomposing the input network into a set of quartets
    topology=get_quartets(net,p)
    
    #---starting values---#
    #starting value for theta
    starting_theta=get_start_theta(topology,
                                    lower_bound=lower_bound,
                                    factor=factor,
                                    tolerance=tolerance)
    if starting_theta<1e-5 
        @error("An error occurred during determining the start value of theta. Setting it to 1e-5 arbitrarily.")
        starting_theta=1e-5 
    end

    #starting values for gamma
    gammas_in_net, nodes_of_hybrid_edge_in_net=get_start_gamma_info(net)
    how_many_gammas=2*net.numHybrids
    length(gammas_in_net)==length(nodes_of_hybrid_edge_in_net) &&
    net.numHybrids==Int(length(gammas_in_net)/2) &&
    how_many_gammas==length(gammas_in_net) || @error("An error occurred during determining the start value(s) of gamma.")
    
    #starting values for tau
    average_momest=get_average_moment_branch_length(topology)
    #tobe fixed: for some reason average_momest sometimes contains one or more negative values, 
    #               but it shouldn't; although it will not affect the final result.
    #all(>=(0), average_momest) || @error("An error occurred during determining the start value(s) of tau. Strating theta seems infeasible.")
    #println(average_momest)
    all_pc_pairs=get_parent_child(net)
    length(all_pc_pairs)==2*net.numHybrids || @error("An error occurred during parametrizing edge lengths.")
    
    #---parameter transformations---#
    #each parameter set contains transformed parameters in the order of [taus,gamma,theta].
    for parent_child in all_pc_pairs
        #params=Float64[] #initializing a parameter set
        #for taus
        params=get_initial_taus(average_momest,parent_child)
        #for gammas
        if isTree(net)
            push!(params,asin(sqrt(1.0)))    
        else    
            for i in 1:how_many_gammas
                push!(params,asin(sqrt(gammas_in_net[i])))    
                i+=2
            end
        end
        #for theta
        push!(params,log(starting_theta))
        #push everything
        push!(parameter_sets,params)
    end

    #---the obejctive function---#
        function objective(params::Array)
            lParam=length(params)
            #backtransform the parameters 
            #taus
            taus=backTransform(average_momest,params,all_pc_pairs[which_pair])
            #gammas
            gammas=Float64[]
            if isTree(net) push!(gammas,1)
            else
                for i in (lParam-net.numHybrids):(lParam-1)
                    push!(gammas,sin(params[i])^2)
                    push!(gammas,1-sin(params[i])^2)
                end
            end
            #theta
            theta=exp(params[lParam])             

            #probabilities for the parental trees
            prob_for_parental_trees=Float64[]
            if net.numHybrids<=1
                prob_for_parental_trees=gammas 
            else
                for i in 1:2
                    for j in 3:length(gammas)
                        push!(prob_for_parental_trees,(gammas[i]*gammas[j]))
                        #this is based on the observation but might not always be true
                    end
                end
            end
                        
            #insert the transformed parameters into the Network object for CL optimization
            topology.theta=theta
            for each_quartet in topology.quartet
                each_quartet.gamma=prob_for_parental_trees[(each_quartet.displayed_tree)]
                each_quartet.average_mom_est_bl=(taus[each_quartet.ntau[1]],taus[each_quartet.ntau[2]],taus[each_quartet.ntau[3]])
            end
            
            #compute negative log composite likelihood using the updated Network object
            neg_log_composite_likelihood=get_negative_log_clikelihood(topology)#Then compute pseudolikelihood
            
            return neg_log_composite_likelihood            
        end

    #---optimization---#
    #some preamble for optimization
    which_pair=0
    res=nothing
    cLikelihood=Inf
    net_after_update=deepcopy(net)

    for params in parameter_sets
        which_pair+=1
        res=Optim.optimize(objective,params,BFGS(linesearch=LineSearches.BackTracking()),Optim.Options(iterations=number_of_itera)) 
        if !isnan(res.minimum) && res.minimum<cLikelihood 
            cLikelihood=res.minimum
            taus,gammas,nodes_of_hybrid_edge_in_net=backtransform_parameters(average_momest,res.minimizer,all_pc_pairs[which_pair],net.numHybrids,nodes_of_hybrid_edge_in_net)
            net_after_update=update_topology(net,taus,gammas,nodes_of_hybrid_edge_in_net)
        end           
    end
    
    return res, net_after_update
end

"""
    update_topology(net_before_update::HybridNetwork,Taus, theta, gammas, nodes_of_hybrid_edge_in_net)

Update optimized parameters on topology
"""
function update_topology(net_before_update::HybridNetwork,
                        taus::Vector{Float64}, 
                        gammas::Vector{Float64}, 
                        nodes_of_hybrid_edge_in_net::Vector{Tuple{Int64, Int64}})
        
    net=deepcopy(net_before_update)
    hnodes=net.hybrid

    for e in net.edge
        u=GetParent(e)
        v=GetChild(e)
        if !(u.hybrid) && !(v.hybrid) # neither u nor v in e=(u,v) are hybrid nodes
            if !(u.hybrid) && v.leaf
                e.length=taus[tauNum(u.number)]
            elseif !(u.hybrid) && !(v.hybrid)
                e.length=taus[tauNum(u.number)]-taus[tauNum(v.number)] 
            end    
        else #at least one of u nor v in e=(u,v) are hybrid nodes
            #only update gamma; reticulation branch lengths are updated next
            for i in 1:length(nodes_of_hybrid_edge_in_net)
                if nodes_of_hybrid_edge_in_net[i]==(u.number,v.number)
                    e.gamma=gammas[i]
                end
            end
        end
    end
    
    for hn in hnodes
        length(hn.edge)==3 || @error("Hybrid node $(hn.number) is not a degree-3 node.")
        
        hnage=0.0
        if !(hn.edge[1].hybrid)
            p1e=taus[abs(GetParent(hn.edge[2]).number)-1]
            p2e=taus[abs(GetParent(hn.edge[3]).number)-1]
            if p1e<p2e 
                (hn.edge[2]).length=0.0
                (hn.edge[3]).length=p2e-p1e
                hnage=p1e
            elseif p1e>p2e 
                (hn.edge[2]).length=p1e-p2e
                (hn.edge[3]).length=0.0
                hnage=p2e
            else
                (hn.edge[2]).length=(hn.edge[3]).length=0.0
                hnage=p1e
            end

            child=GetChild(hn.edge[1])
            if child.leaf
                hn.edge[1].length=hnage
            elseif child.hybrid
                @error "Descendant of a hybrid node cannot be another hybrid node."
            else
                hn.edge[1].length=hnage-taus[abs(child.number)-1]
            end
            
            
        elseif !(hn.edge[2].hybrid)
            p1e=taus[abs(GetParent(hn.edge[1]).number)-1]
            p2e=taus[abs(GetParent(hn.edge[3]).number)-1]
            if p1e<p2e 
                (hn.edge[1]).length=0.0
                (hn.edge[3]).length=p2e-p1e
                hnage=p1e
            elseif p1e>p2e 
                (hn.edge[1]).length=p1e-p2e
                (hn.edge[3]).length=0.0
                hnage=p2e
            else
                (hn.edge[1]).length=(hn.edge[3]).length=0.0
                hnage=p1e
            end

            child=GetChild(hn.edge[2])
            if child.leaf
                hn.edge[2].length=hnage
            elseif child.hybrid
                @error "Descendant of a hybrid node cannot be another hybrid node."
            else
                hn.edge[2].length=hnage-taus[abs(child.number)-1]
            end
        else 
            p1e=taus[abs(GetParent(hn.edge[1]).number)-1]
            p2e=taus[abs(GetParent(hn.edge[2]).number)-1]
            if p1e<p2e 
                (hn.edge[1]).length=0.0
                (hn.edge[2]).length=p2e-p1e
                hnage=p1e
            elseif p1e>p2e 
                (hn.edge[1]).length=p1e-p2e
                (hn.edge[2]).length=0.0
                hnage=p2e
            else
                (hn.edge[1]).length=(hn.edge[2]).length=0.0
                hnage=p1e
            end

            child=GetChild(hn.edge[3])
            if child.leaf
                hn.edge[3].length=hnage
            elseif child.hybrid
                @error "Descendant of a hybrid node cannot be another hybrid node."
            else
                hn.edge[3].length=hnage-taus[abs(child.number)-1]
            end
        end


    end

    #set the number of digits at mantissa and round
    #numdigits is a global variable and set as 5 by default (in Objects.jl)
    for e in net.edge
        e.length=round(e.length,digits=numdigits)
        if e.gamma<1 e.gamma=round(e.gamma,digits=numdigits) end
    end
    
    return net
end



"""
    get_start_gamma_info(net::HybridNetwork)

Returns the inheritance probabilities in the network topology\\
and the parent and child nodes of the branches the inheritance probability is assigned.
"""
function get_start_gamma_info(net::HybridNetwork)
    gammas_in_net=Float64[]
    nodes_of_hybrid_edge_in_net=Tuple{Int64, Int64}[]
    
    for each_edge in net.edge
        if (each_edge.hybrid)
            push!(gammas_in_net,each_edge.gamma)
            head_and_tail=Int64[]
            for the_two_nodes in each_edge.node
                push!(head_and_tail,the_two_nodes.number)
            end
            #head_and_tail should contain two node numbers, in the order of head and tail.
            push!(nodes_of_hybrid_edge_in_net,(head_and_tail[2],head_and_tail[1]))
        end
    end
    
    return gammas_in_net, nodes_of_hybrid_edge_in_net
end




"""
    get_parent_child

Function to get a vector of the parent-child relationships of node numbers using the node numbers assigned 
in HybridNetwork then transform it into PhyNEST's numbering system (e.g., node number -2 represents the root in 
HybridNetwork but we use 1). The length of the vector must be equal to the number of speciation times in the
topology (or the number of non-leaf vertices) that equals to J=n+h-1 where n represents the number of leaves
and h represents the number of reticulations.

### Input
**`net`**       A tree/network in Type object PhyloNetworks.HybridNetwork\\
"""
function get_parent_child(net::HybridNetwork) 
    parentchildss=Array[]
    dTrees=displayedTrees(net,0.0,nofuse=true)
    #for dt in dTrees
    #    display(PhyloNetworks.printEverything(dt))
    #end

    for dt in dTrees
        pc=[]
        for n in dt.node
            i=n.number
            while i<0
                parent=parentnode(n.number,dt)
                if parent===nothing #the node is root
                    push!(pc,[0,i])
                elseif parent>0
                    parent=parentnode(parent,dt)
                    push!(pc,[parent,i])
                else
                    push!(pc,[parent,i])
                end
                break
            end
        end
        push!(parentchildss,pc)
    end

#=
    for dt in dTrees        
        for n in dt.node
            i=n.number
            while i<0 #the node is not a leaf
                parent=parentnode(n.number,dt)
                if parent===nothing #the node is root
                    push!(parentchild,[0,i])
                else
                    push!(parentchild,[parent,i])
                end
                break
            end
        end
    end
=#
    for parentchild in parentchildss    
        parentchild=unique(parentchild)
        #organize array, root on top and moves down the network
        parentchild=sort(parentchild, rev=true)
    end
    
    for parentchild in parentchildss
        for element in parentchild #just renaming for PhyNe from the number assigned from PhyloNetworks
            i=element[1] #parent
            j=element[2] #child
            if i!==0 
                element[1]=abs(i)-1
            else
                element[1]=0
            end
            #basically j is always non-zero
            #if j!==0
                element[2]=abs(j)-1
            #else
            #    element[2]=0
            #end
        end
    end

    npcs=[]
    for pc in parentchildss
        npc=sort(pc, rev=false)
        push!(npcs,npc)
    end
    parentchildss=npcs
    
    return parentchildss
end

"""
    get_initial_taus(ave::Array, parentchild::Array)

Returns the ratio node ages that are going to be used for the unconstrained optimization.\\
A few arbitrary assumptions were made in case unrealistic values occurs... probably okay.
"""
function get_initial_taus(ave::Array, parentchild::Array)
    
    transTau=Float64[]

    for element in parentchild
        parent=element[1]
        child=element[2]
        if parent==0
            try
                push!(transTau,log(ave[child]))
            catch e 
                push!(transTau,log(20))
            end
        else
            if ave[child]/ave[parent] > 1
                push!(transTau,asin(sqrt(0.999)))
            elseif ave[child]/ave[parent] == 1
                push!(transTau,asin(sqrt(0.999)))
            else
                try
                    push!(transTau,asin(sqrt(ave[child]/ave[parent])))
                catch e push!(transTau,asin(sqrt(0.999)))
                end
            end
        end
    end

    return transTau

end

"""
    backTransform(ave::Array, params::Array,pc::Array)

Back transforms the transformed node ages with a help of couple more information...
"""
function backTransform(ave::Array, params::Array,pc::Array)
    TausMat=zeros(Float64, 1, length(ave))
    Taus=Float64[]
    pcc=0
    for element in pc
        pcc+=1
        #println("ElementParentChild: $element")
        parent=element[1]
        child=element[2]
        if parent===0
            TausMat[child]=exp(params[child])
            #append!(Taus,exp(params[child]))
        else
            TausMat[child]=TausMat[parent]*(sin(params[pcc])^2)
            #append!(Taus,(Taus[parent]*(sin(params[child])^2)))
        end
    end
    
    for n in 1:size(TausMat)[2]
        append!(Taus,TausMat[n])
    end

    return Taus
end

"""
    backtransform_parameters(average_momest::Array, res_minimizer::Array, pc::Array,how_many_hybrid_nodes::Int64,nodes_of_hybrid_edge_in_net::Vector{Tuple{Int64, Int64}})

Currently just back transform the optimized parameter set into interpretable values, \\
soon this function will update the input topology parameters using these values.
"""
function backtransform_parameters(average_momest::Array, 
                                    res_minimizer::Array, pc::Array,
                                    how_many_hybrid_nodes::Int64,
                                    nodes_of_hybrid_edge_in_net::Vector{Tuple{Int64, Int64}})

    Taus=backTransform(average_momest,res_minimizer,pc)
    #theta=exp(res_minimizer[length(res_minimizer)])*2
    gammas=Float64[]
        if how_many_hybrid_nodes==0
            push!(gammas,1)
        else# how_many_hybrid_nodes>0
            for i in (length(res_minimizer)-(how_many_hybrid_nodes)):(length(res_minimizer)-1)
                push!(gammas,sin(res_minimizer[i])^2)
                push!(gammas,1-sin(res_minimizer[i])^2)
            end
        end
    #println("gammas=$(gammas)")
    #println("nodes_of_hybrid_edge_in_net=$(nodes_of_hybrid_edge_in_net)")
    #println("Taus=$(Taus)")
    #println("theta=$(exp(res_minimizer[length(res_minimizer)])*2)")
    
    return Taus, gammas, nodes_of_hybrid_edge_in_net
end

