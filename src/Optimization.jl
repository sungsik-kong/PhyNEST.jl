#written by Sungsik Kong 2021-2022
"""
    do_optimization(net::HybridNetwork, p::Phylip;
                        lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64,
                        number_of_itera=1000::Int64)

Optimizes the composite likelihood of the network topology given the data.
"""
function do_optimization(net::HybridNetwork, p::Phylip;
                        lower_bound=0.00001::Float64,
                        factor=2.0::Float64,
                        tolerance=0.01::Float64,
                        number_of_itera=1000::Int64,
                        update_parameters=true::Bool)
    
    params=Float64[]
    N=get_quartets(net,p)
    new_topology=N
    
    #get the starting values 
    theta=get_start_theta(N,lower_bound=lower_bound,factor=factor,tolerance=tolerance)
    
    average_momest=get_average_moment_branch_length(N)
    pc=get_parent_child(net)
    
    gammas_in_net, nodes_of_hybrid_edge_in_net=get_start_gamma_info(net)
    how_many_hybrid_nodes=net.numHybrids
    net.numHybrids==length((gammas_in_net)/2)
    how_many_gammas=length(gammas_in_net)
    
    #trasform the starting parameters and push them to the params array
    initTau=get_initial_taus(average_momest,pc)
    for element in initTau
        push!(params,element)
    end

    if gammas_in_net==1
        push!(params,asin(sqrt(1.0)))    
    else    
        for i in 1:how_many_gammas
            push!(params,asin(sqrt(gammas_in_net[i])))    
            i+=2
        end
    end
  
    push!(params,log(theta))
    
    #the obejctive function
    function objective(params::Array)
        #backtransform the parameters 
        theta=exp(params[length(params)])  
        if theta<1e-5 theta=1e-5 end

        Taus=backTransform(average_momest,params,pc)

        gammas=Float64[]
        if how_many_hybrid_nodes==0
            push!(gammas,1)
        else
            for i in (length(params)-(how_many_hybrid_nodes)):(length(params)-1)
                push!(gammas,sin(params[i])^2)
                push!(gammas,1-sin(params[i])^2)
            end
        end

        prob_for_parental_trees=Float64[]
        if how_many_hybrid_nodes<=1
            prob_for_parental_trees=gammas 
        else
            for i in 1:2
                for j in 3:length(gammas)
                    push!(prob_for_parental_trees,(gammas[i]*gammas[j]))
                    #this is based on the observation but might not always be true
                end
            end
        end
        
        #reimplement the transformed parameters into the Network object for CL optimization
        new_topology.theta=theta

        for each_quartet in new_topology.quartet
            each_quartet.gamma=prob_for_parental_trees[(each_quartet.displayed_tree)]
        end

        for each_quartet in new_topology.quartet
            each_quartet.average_mom_est_bl=(Taus[each_quartet.ntau[1]],Taus[each_quartet.ntau[2]],Taus[each_quartet.ntau[3]])
        end
        
        #compute negative log composite likelihood using the updated Network object
        neg_log_composite_likelihood=get_negative_log_clikelihood(new_topology)#Then compute pseudolikelihood
        
        return neg_log_composite_likelihood
        
    end

    #optimize the objective function and kind of visualize the minimizer in interpretable values
    res=Optim.optimize(objective,params,BFGS(linesearch=LineSearches.BackTracking()),Optim.Options(iterations = number_of_itera))
    net_after_update=deepcopy(net)
    if(update_parameters)
        Taus, theta, gammas, nodes_of_hybrid_edge_in_net=backtransform_parameters(average_momest, res.minimizer, pc,how_many_hybrid_nodes,nodes_of_hybrid_edge_in_net)
        net_after_update=update_topology(net,Taus, theta, gammas, nodes_of_hybrid_edge_in_net)
        #println(Taus, theta, gammas, nodes_of_hybrid_edge_in_net)
    end
    
    return res, net_after_update
end

"""
    get_start_gamma_info(net::HybridNetwork)

Returns the inheritance probabilities in the network topology\\
and the parent and child nodes of the branches the inheritance probability is assigned.
"""
function get_start_gamma_info(net::HybridNetwork)
    gammas_in_net=Float64[]
    nodes_of_hybrid_edge_in_net=Tuple{Int64, Int64}[]
    edges=net.edge

    for each_edge in edges
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
in HybridNetwork then transform it into our numbering system (e.g., node number -2 represents the root in 
HybridNetwork but we use 1). The length of the vector must be equal to the number of speciation times in the
topology (or the number of non-leaf vertices) that equals to J=n+h-1 where n represents the number of leaves
and h represents the number of reticulations.

### Input
**`net`**       A tree/network in Type object PhyloNetworks.HybridNetwork\\
"""
function get_parent_child(net::HybridNetwork) 
    parentchild=Array[]
    for n in net.node
        i=n.number
        while i<0 #the node is not a leaf
            parent=parentnode(n.number,net)
            if parent===nothing #the node is root
                push!(parentchild,[0,i])
            else
                push!(parentchild,[parent,i])
            end
            break
        end
    end
 
    #organize array, root on top and moves down the network
    parentchild=sort(parentchild, rev=true)
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

    return parentchild
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
    theta=exp(res_minimizer[length(res_minimizer)])*2
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
    
    return Taus, theta, gammas, nodes_of_hybrid_edge_in_net
end

"""
    update_topology(net_before_update::HybridNetwork,Taus, theta, gammas, nodes_of_hybrid_edge_in_net)

Update optimized parameters on topology
"""
function update_topology(net_before_update::HybridNetwork,Taus, theta, gammas, nodes_of_hybrid_edge_in_net)
    #printEdges(net)
    #println(writeTopologyLevel1(net))
    #println(Taus, theta, gammas, nodes_of_hybrid_edge_in_net)
    #println(Taus)
    net=deepcopy(net_before_update)
    n_taus=length(Taus)
    
    for branch in net.edge
        parent=GetParent(branch)
        child=GetChild(branch)
        
        parent_node=parent.number
        child_node=child.number
        
        tau_num_parent=tauNum(parent_node)
        tau_num_child=tauNum(child_node)
        
        #branch lengths
        if tau_num_parent<=n_taus #&& tau_num_child<=n_taus
            if child_node<0 && !(branch.hybrid)   # tree branch
                branch.length=Taus[tau_num_parent]-Taus[tau_num_child]
            elseif child_node>0 && !(branch.hybrid)   #tree branch leading to the tip
                #branch lengths for the hybrid edge are set as 0.0
                branch.length=Taus[tau_num_parent]
            elseif child_node>0 && (branch.hybrid)   #hybrid branch leading to the hybrid node
                branch.length=Taus[tau_num_parent]
            else #child_node<0 && (branch.hybrid)   
                branch.length=0.0 
            end
        else
            branch.length=10^-10
        end

        #gamma
        for i in 1:length(nodes_of_hybrid_edge_in_net)
            if (parent_node,child_node)==nodes_of_hybrid_edge_in_net[i]
                this_gamma=gammas[i]
                if this_gamma==NaN
                    this_gamma=1/length(nodes_of_hybrid_edge_in_net)
                    branch.gamma=this_gamma
                else
                    branch.gamma=this_gamma
                end
            elseif (child_node,parent_node)==nodes_of_hybrid_edge_in_net[i]
                this_gamma=gammas[i]
                    if this_gamma==NaN
                    this_gamma=1/length(nodes_of_hybrid_edge_in_net)
                    branch.gamma=this_gamma
                else
                    branch.gamma=this_gamma
                end
            else 
                continue
            end
        end

    end
    
    #fix hybrid branch lengths
    hybrid_nodes=net.hybrid
    for each_hybrid_node in hybrid_nodes
        edge1=each_hybrid_node.edge[1]
        edge2=each_hybrid_node.edge[2]
        edge3=each_hybrid_node.edge[3]

        edge1_length=edge1.length
        edge2_length=edge2.length
        edge3_length=edge3.length

        the_three_edge_lengths=[edge1_length,edge2_length,edge3_length]
        #println([edge1_length,edge2_length,edge3_length])
        sorted=sort([edge1_length,edge2_length,edge3_length])
        shortest_length=minimum(the_three_edge_lengths)

        if edge1.hybrid==false
            if shortest_length!==0
                edge1.length=shortest_length
            else
                edge1.length=sorted[2]
            end
        elseif edge1.hybrid==true
            if edge1.length==shortest_length
                edge1.length=0
            else
                edge1.length=edge1_length-shortest_length
            end
        end

        if edge2.hybrid==false
            edge2.length=shortest_length
        elseif edge2.hybrid==true
            if edge2.length==shortest_length
                edge2.length=0
            else
                edge2.length=edge2_length-shortest_length
            end
        end

        if edge3.hybrid==false
            edge3.length=shortest_length
        elseif edge3.hybrid==true
            if edge3.length==shortest_length
                edge3.length=0
            else
                edge3.length=edge3_length-shortest_length
            end
        end        
    end

    return net
end