global const numdigits = 3
#written by Sungsik Kong 2021-2022
"""
    do_optimization(net::HybridNetwork, p::Phylip;
                    lower_bound=0.00001::Float64,
                    factor=2.0::Float64,
                    tolerance=0.01::Float64,
                    number_of_itera=1000::Int64,
                    minimum_theta=1e-5::Float64)

Optimizes the composite likelihood of the network topology given the data.


"""
function do_optimization(initial_topology::HybridNetwork, p::Phylip;
                        lower_bound=0.00001::Float64,
                        factor=2.0::Float64,
                        tolerance=0.01::Float64,
                        number_of_itera=1000::Int64,
                        minimum_theta=1e-5::Float64)
    
    #---check topology---#
    #weather or not gamma values are specified, if not, arbitrarily insert 0.5.
    for e in initial_topology.edge
        if e.hybrid && !(0<e.gamma<1) e.gamma=0.5 end
    end
    #check branch lengths to if I want to use it
    #maany other stuffs can go in here if needed
    
    #---preamble---#
    #topology of interest
    #we use deepcopy so it does not mess up the original topology
    net=deepcopy(initial_topology)
    #decomposing the input network into a set of quartets
    topology=get_quartets(net,p)
    #length should be 2*h because we parametrize node ages two ways as we cannot identify reticulation node age
    parameter_sets=Vector{Array{Float64}}()
    
    #---starting values---#
    #starting value for theta
    starting_theta=get_start_theta(topology, lower_bound=lower_bound, factor=factor, tolerance=tolerance)
    if starting_theta<minimum_theta starting_theta=minimum_theta end #estimated theta should be great than setted minimum theta
    #starting values for tau
    average_momest=get_average_moment_branch_length(topology)
    #starting values for gamma
    gammas_in_net, nodes_of_hybrid_edge_in_net=get_start_gamma_info(net)
    #get all parent-child relationships
    all_pc_pairs=get_parent_child(net)
    
    #---parameter transformations---#
    #each parameter set contains transformed parameters in the order of [taus,gamma,theta].
    for parent_child in all_pc_pairs
        #params=Float64[] #initializing a parameter set
        #for taus
        params=transform_taus(average_momest,parent_child)
        #for gammas
        if isTree(net)
            push!(params,asin(sqrt(1.0)))    
        else    
            for i in 1:(2*net.numHybrids)
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
            lengthavemom=length(average_momest)
            #backtransform the parameters 

            #taus
            taus=backtransform_taus(lengthavemom,params,all_pc_pairs[which_pair])
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
    how_many_hybrid_nodes=0
    res=nothing
    cLikelihood=Inf
    net_after_update=deepcopy(net)
    gammas=Float64[]

    for params in parameter_sets
        which_pair+=1
        res=Optim.optimize(objective,params,BFGS(linesearch=LineSearches.BackTracking()),Optim.Options(iterations=number_of_itera)) 
        if !isnan(res.minimum) && res.minimum<cLikelihood 
            cLikelihood=res.minimum
            lengthavemom=length(average_momest)
            taus=backtransform_taus(lengthavemom,res.minimizer,all_pc_pairs[which_pair])
            
            if isTree(net) push!(gammas,1)
            else# how_many_hybrid_nodes>0
                for i in (length(res.minimizer)-(how_many_hybrid_nodes)):(length(res.minimizer)-1)
                    push!(gammas,sin(res.minimizer[i])^2)
                    push!(gammas,1-sin(res.minimizer[i])^2)
                end
            end
            
            #taus,gammas,nodes_of_hybrid_edge_in_net=backtransform_parameters(average_momest,res.minimizer,all_pc_pairs[which_pair],net.numHybrids,nodes_of_hybrid_edge_in_net)            
            net_after_update=update_topology(net,taus,gammas,nodes_of_hybrid_edge_in_net)
        end           
    end
    
    return res, net_after_update
end


function isleaf(n::Node)
    if n.leaf return true
    else return false
    end
end

"""
    update_topology(net_before_update::HybridNetwork,Taus, theta, gammas, nodes_of_hybrid_edge_in_net)

Update optimized parameters on topology
"""
function update_topology(net_before_update::HybridNetwork,
                        taus::Vector{Float64}, 
                        gammas::Vector{Float64}, 
                        nodes_of_hybrid_edge_in_net::Vector{Tuple{Int64, Int64}})
        
    #---Preamble---#
    net=deepcopy(net_before_update)
    
    #---#initialization---#
    for e in net.edge
        e.length=-1.0
        e.gamma=-1.0
    end

    #--tree node ages and gamma---#
    for e in net.edge
        u=GetParent(e) #parent node
        v=GetChild(e) #child node
        if !(u.hybrid) && !(v.hybrid) # neither u nor v in e=(u,v) are hybrid nodes
            if !(u.hybrid) && v.leaf #treminal tree branch that does not starts at a hybrid node
                e.length=taus[tauNum(u.number)]
            else#if !(u.hybrid) && !(v.hybrid)
                e.length=taus[tauNum(u.number)]-taus[tauNum(v.number)] 
            end    
        elseif !(u.hybrid) && (v.hybrid)
            #only update gamma; reticulation branch lengths are updated soon
            for i in 1:length(nodes_of_hybrid_edge_in_net)
                if nodes_of_hybrid_edge_in_net[i]==(tauNum(u.number),tauNum(v.number))
                    e.gamma=gammas[i]
                end
            end
        end
    end
    
    #---reticulation node age---#
    for n in net.node
        if n.hybrid
            nodeage=Inf
            #get the node age of the hybrid node here
            for e in n.edge
                u=GetParent(e)
                v=GetChild(e)
                if (v.hybrid)
                    if taus[tauNum(u.number)]<nodeage
                        nodeage=taus[tauNum(u.number)]
                    else continue
                    end
                else
                    continue
                end                
            end

            #assigning length to the edges involving reticulation edge
            for e in n.edge
                u=GetParent(e)
                v=GetChild(e)
                if(isleaf(v)) #child is leaf
                    e.length=nodeage
                elseif v.hybrid #reticulation edge
                    e.length=taus[tauNum(u.number)]-nodeage
                else #child node is a tree node
                    e.length=nodeage-taus[tauNum(v.number)]
                end
            end
        end
    end

    #---Postprocessing---#
    for e in net.edge
        if e.length<0 
            e.length=0.0 #setting branch length to zero when it is negative
        else 
            e.length=round(e.length,digits=numdigits) 
        end

        if e.gamma<1 e.gamma=round(e.gamma,digits=numdigits) end
    end
    
    return net
end


"""
    get_initial_taus(ave::Array, parentchild::Array)

Returns the ratio node ages that are going to be used for the unconstrained optimization.\\
A few arbitrary assumptions were made in case unrealistic values occurs... probably okay.
"""
function transform_taus(ave::Array, parentchild::Array; missingrootage=3.0::Float64, nearlyone=0.999::Float64)
    
    transformedtau=Vector{Float64}()

    for pair in parentchild
        u=pair[1]
        v=pair[2]
        if u==0 #root
            try
                push!(transformedtau,log(ave[v]))
            catch 
                push!(transformedtau,log(maximum(ave)+missingrootage))
            end
        else #nonroot
            if ave[v]/ave[u] >= 1
                push!(transformedtau,asin(sqrt(nearlyone)))
            else
                try
                    push!(transformedtau,asin(sqrt(ave[v]/ave[u])))
                catch 
                    push!(transformedtau,asin(sqrt(nearlyone)))
                end
            end
        end
    end

    return transformedtau
end

"""
    backTransform(ave::Array, params::Array,pc::Array)

Back transforms the transformed node ages with a help of couple more information...
"""
function backtransform_taus(lengthavemom::Integer,params::Array,pc::Array)
    pcc=0
    taus=zeros(lengthavemom)

    for pair in pc
        pcc+=1
        u=pair[1]
        v=pair[2]

        if u==0
            taus[v]=exp(params[v])
        else
            taus[v]=taus[u]*(sin(params[pcc])^2)
        end
    end

    return taus
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
    displayedtrees=displayedTrees(net, 0.0, nofuse=true)

    for dt in displayedtrees
        pc=[]
        for n in dt.node
            i=n.number
            while i<0 #while the node is a tree node
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
        
        for element in pc #just renaming for PhyNe from the number assigned from PhyloNetworks
            i=element[1] #parent
            j=element[2] #child
            if i!==0 
                element[1]=tauNum(i)#abs(i)-1
            else
                element[1]=0
            end
            #basically j is always non-zero
            #if j!==0
                element[2]=tauNum(j)#abs(j)-1
            #else
            #    element[2]=0
            #end
        end

        pc=unique(pc)
        pc=sort(pc, rev=false)

        push!(parentchildss,pc)
    end
    
    return parentchildss
end















"""
    get_start_gamma_info(net::HybridNetwork)

Returns the inheritance probabilities in the network topology\\
and the parent and child nodes of the branches the inheritance probability is assigned.
"""
function get_start_gamma_info(net::HybridNetwork)
    gammas_in_net=Float64[]
    nodes_of_hybrid_edge_in_net=Tuple{Int64, Int64}[]
    
    for e in net.edge
        if(e.hybrid)
            #head_and_tail=Int64[]
            u=GetParent(e).number #parent node
            v=GetChild(e).number #child node
            the_two_nodes=(tauNum(u),tauNum(v))
            #for the_two_nodes in each_edge.node
            #    println(the_two_nodes)
            #    println("---")
            #push!(head_and_tail,the_two_nodes.number)
            #end
            push!(gammas_in_net,e.gamma)
            #head_and_tail should contain two node numbers, in the order of head and tail.
            push!(nodes_of_hybrid_edge_in_net,the_two_nodes)
        end
    end
    
    return gammas_in_net, nodes_of_hybrid_edge_in_net
end



