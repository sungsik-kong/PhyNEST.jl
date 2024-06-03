#written by Sungsik Kong 2021-2022
global const alpha=4/3

"""
    get_start_theta(N::Network; lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
Using the Network object that contains information on quartets extracted from the topology and\\
site pattern frequencies from the data, estimates the `reasonable' theta value that can be\\
used as the starting point for the composite likelihood optimization. \\
Approximate interval where the true theta may lie is estimated through two procedures:\\
somewhat loose interval that keeps the branch lengths positive from the moment estimator,\\
and then further tightens the interval using the golden section seach.
"""
function get_start_theta(net::Network; lower_bound=1e-5::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
    
    function find_theta(theta::Float64)
        new_topology=deepcopy(net)
        new_topology.theta=theta

        for each_quartet in new_topology.quartet
            each_quartet.momestlength=methodofmomentestimator(each_quartet,theta)
        end
        average_momest=get_average_moment_branch_length(new_topology)
        #println(average_momest)
        for each_quartet in new_topology.quartet
            each_quartet.average_mom_est_bl=(average_momest[each_quartet.ntau[1]],average_momest[each_quartet.ntau[2]],average_momest[each_quartet.ntau[3]])
        end
        
        neg_log_composite_likelihood=get_negative_log_clikelihood(new_topology)#Then compute pseudolikelihood
        
        return neg_log_composite_likelihood
    end

    lower,upper=get_upper_lower_theta(net, lower_bound=lower_bound, factor=factor)
    new_lower,new_upper=golden_section_search(find_theta,lower,upper,tolerance)
    res=Optim.optimize(find_theta,new_lower,new_upper)
    theta=res.minimizer

    return theta
end

function get_start_theta(net::HybridNetwork, p::Phylip; lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
    setofquartets=get_quartets(net,p)
    theta=get_start_theta(setofquartets; lower_bound=lower_bound,factor=factor,tolerance=tolerance)
    return theta
end

"""
    get_negative_log_clikelihood(N::Network)

Computes the negative log coposite likelihood from the Network object. It first computes\\
the true site pattern probabilities for each quartet extracted from the network,\\
which can be computed as sum of all the true site pattern probabilities for a quartet extracted from\\
each parental tree multiplies by the inheritance probability assigned to that quartet (i.e., the\\
probability that is assigned to the parental tree that the quartet is extracted from.\\
Then, it computes quartet likelihood, followed by multiplification of the likelihoods\\
for every quartet that is extracted from the network. This will be a very small number, so we\\
use BigFloat(). Finally we cat Float64 value by multiplying -1.0 for the log of that value. 
"""
function get_negative_log_clikelihood(N::Network)
    
    theta=N.theta

    #get prob for each quartet
    probs=[]
    how_many_trees=length(N.gamma)
    all_quartets=N.quartet
    how_many_quartets=Int(length(N.quartet)/how_many_trees)

    #---get true site pattern probabilities for each quartet---#
    for i in 1:how_many_quartets
        site_pattern_prababilities=zeros(15)#[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        this_quartet=all_quartets[i].quartet
        for each_quartet in all_quartets
            if each_quartet.quartet==this_quartet
                site_pattern_prababilities+=get_true_probabilities(each_quartet,theta,alpha)
            else continue
            end
        end
        site_pattern_prababilities = (weights .* site_pattern_prababilities)
        push!(probs,site_pattern_prababilities)
    end

    #---computing quartet likelihood---#
    quartetlikelihoods=[]
    for i in 1:how_many_quartets
        site_pattern_frequencies=all_quartets[i].mspcountsNET
        site_pattern_probs=probs[i]
        if (sum(ismissing.(site_pattern_probs))!=0) 
            #likelihood=-10^308
            likelihood=Inf
        elseif (sum(any(t->t<=0,site_pattern_probs))==0)
            #likelihood=BigFloat(1.0)
            likelihood=0.0
            for i in 1:15
                #likelihood*=BigFloat(site_pattern_probs[i])^site_pattern_frequencies[i]
                likelihood+=-1*site_pattern_frequencies[i]*log(site_pattern_probs[i])
            end
        else 
            #likelihood=-10^308
            likelihood=Inf
        end 
        #push!(qliks,likelihood)
        push!(quartetlikelihoods,likelihood)
    end

    #get clikelihood
    #clikelihood=prod(qliks)
    clikelihood=sum(quartetlikelihoods)
    #make it negative log and back to Float from BigFloat
    #neg_log_clik=Float64(-1.0*log(clikelihood))

    upd_log_clik=sum(clikelihood)
    #println(clikelihood)

    #println("------")
    #println("neg_log_clik: $neg_log_clik")
    #println(neg_log_clik==upd_log_clik)
    #println("------")
    #if isnan(upd_log_clik) println("$upd_log_clik is nan!") 
    #    upd_log_clik=Inf
    #end
    #println("$upd_log_clik")

    return upd_log_clik
end

"""
    get_true_probabilities(each_quartet::quartets,theta::Float64,alpha::Float64)

Computes the tru site pattern probabilities for each quartet in a network,\\
multiplied by the inheritance probability of that quartet.
"""
function get_true_probabilities(each_quartet::quartets,theta::Float64,alpha::Float64)
    
    t1=each_quartet.average_mom_est_bl[1]
    t2=each_quartet.average_mom_est_bl[2]
    t3=each_quartet.average_mom_est_bl[3]

    if each_quartet.symtype==0 
        probabilities=each_quartet.gamma*GetTrueProbsSymm(t1,t2,t3,theta,alpha)
    else probabilities=each_quartet.gamma*GetTrueProbsAsymmTypes(each_quartet.symtype,t1,t2,t3,theta,alpha)
    end
    
    return probabilities
end


"""
    get_upper_lower_theta(N::Network;lower_bound=0.00001::Float64, factor=2.0::Float64)

Get somewhat loose intherval of theta. In particular, the interval is initially determined as\\
the user specified lower bound (default=0.00001) and the lower bound * user specified factor (dafault=2.0).\\
If all branch lengths estimated using the theta=upper bound happens to be positive, the upper bound\\
is `feasible` theta value, and we iterate this process for multiple iterations (i) until using \\
theta=upper bound becomes infeasible i.e., at least one of the branch lengths estimated for \\
each quartet is negative. The the largest upper bound used is returned.
"""
function get_upper_lower_theta(net::Network;lower_bound=0.00001::Float64, factor=2.0::Float64)
    #---Preamble---#
    iteration=1
    positive=true
    upper_bound=lower_bound*(factor^iteration)

    #---ad hoc---#
    #using theta as the given upper bound, run mom estimates of node ages for all quartets,
    #if all momest node ages are positive, increase theta by twofold;
    #until the theta esitimates any one of the node ages is negative.
    #Once negative, we return the largest theta value tested (i.e., lower*(factor&(interation-1)))
    while(positive)
        for each_quartet in net.quartet
            upper=lower_bound*(factor^iteration)
            t1,t2,t3=methodofmomentestimator(each_quartet,upper)
            if t1>0 && t2>0 && t3>0
                iteration=iteration+1
            else 
                positive=false 
                if iteration==1 upper_bound=0.1
                    break
                else
                    upper_bound=lower_bound*(factor^(iteration-1))
                    break
                end
            end
        end
    end

    return lower_bound,upper_bound
end

"""
    function golden_section_search(f, lower::Float64, upper::Float64, tolerance::Float64)

Using the loose interval for theta estimated from function `get_upper_lower_theta()`, 
golden section search further tighten the interval for the feasible theta using the golden ratio.
"""
function golden_section_search(f, lower::Float64, upper::Float64, tolerance::Float64)
    
    goldenratio=2/(sqrt(5)+1)
    
    #Use the golden ratio to set the initial test points
    val1=upper - goldenratio*(upper-lower)
    val2=lower + goldenratio*(upper-lower)

    #evaluate function at the test points
    func1=f(val1)
    func2=f(val2)

    iteration = 0

    while (abs(upper-lower) > tolerance)
        iteration = iteration+1
        if func2 > func1
            upper=val2
            val2=val1
            func2=func1

            val1=upper-goldenratio*(upper-lower)
            func1=f(val1)
        else
            lower=val1
            val1=val2
            func1=func2

            val2=lower+goldenratio*(upper-lower)
            func2=f(val2)
        end
    end

    return lower,upper
end
