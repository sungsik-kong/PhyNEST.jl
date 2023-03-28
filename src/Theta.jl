#written by Sungsik Kong 2021-2022
#Shame on my bad naming sense...again!
function get_start_theta(net::HybridNetwork, p::Phylip; lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
    N=get_quartets(net,p)
    theta=get_start_theta(N; lower_bound=lower_bound,factor=factor,tolerance=tolerance)
    return theta
end

"""
    get_start_theta(N::Network; 
                    lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
Using the Network object that contains information on quartets extracted from the topology and\\
site pattern frequencies from the data, estimates the `reasonable' theta value that can be\\
used as the starting point for the composite likelihood optimization. \\
Approximate interval where the true theta may lie is estimated through two procedures:\\
somewhat loose interval that keeps the branch lengths positive from the moment estimator,\\
and then further tightens the interval using the golden section seach.
"""
function get_start_theta(N::Network; 
                        lower_bound=0.00001::Float64,factor=2.0::Float64,tolerance=0.01::Float64)
    function find_theta(theta::Float64)
        new_topology=deepcopy(N)
        new_topology.theta=theta

        for each_quartet in new_topology.quartet
            each_quartet.momestlength=momentEstimate(each_quartet,theta)
        end
        average_momest=get_average_moment_branch_length(new_topology)
        for each_quartet in new_topology.quartet
            each_quartet.average_mom_est_bl=(average_momest[each_quartet.ntau[1]],average_momest[each_quartet.ntau[2]],average_momest[each_quartet.ntau[3]])
        end
        
        neg_log_composite_likelihood=get_negative_log_clikelihood(new_topology)#Then compute pseudolikelihood
        
        return neg_log_composite_likelihood
    end

    lower,upper=get_upper_lower_theta(N, lower_bound=lower_bound, factor=factor)
    new_lower,new_upper=golden_section_search(find_theta,lower,upper,tolerance)
    res=Optim.optimize(find_theta,new_lower,new_upper)
    theta=res.minimizer

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
    alpha=4/3

    #get prob for each quartet
    probs=[]
    how_many_trees=length(N.gamma)
    all_quartets=N.quartet
    how_many_quartets=length(N.quartet)/how_many_trees
    how_many_quartets=Int(how_many_quartets)
    for i in 1:how_many_quartets
        site_pattern_prababilities=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        this_quartet=all_quartets[i].quartet
        for each_quartet in all_quartets
            if each_quartet.quartet==this_quartet

                site_pattern_prababilities+=get_true_probabilities(each_quartet,theta,alpha)
            else continue
            end
        end
        weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
        site_pattern_prababilities = (weights .* site_pattern_prababilities)
        push!(probs,site_pattern_prababilities)
    end

    qliks=[]
    for i in 1:how_many_quartets
        site_pattern_frequencies=all_quartets[i].mspcountsNET
        site_pattern_probs=probs[i]
        if (sum(ismissing.(site_pattern_probs))!=0) 
            likelihood=-10^308
        elseif (sum(any(t->t<=0,site_pattern_probs))==0)
            likelihood=BigFloat(1.0)
            for i in 1:15
                likelihood*=BigFloat(site_pattern_probs[i])^site_pattern_frequencies[i]
            end
        else 
            likelihood=-10^308
        end 
        push!(qliks,likelihood)
    end

    #get clikelihood
    clikelihood=prod(qliks)

    #make it negative log and back to Float from BigFloat
    neg_log_clik=Float64(-1.0*log(clikelihood))

    return neg_log_clik
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
function get_upper_lower_theta(N::Network;lower_bound=0.00001::Float64, factor=2.0::Float64)
    i=1
    positive=true
    lower=lower_bound
    upper=lower*((factor)^i)
    
    while (positive=true)
        f=[]
        for each_quartet in N.quartet
            upper=lower*((factor)^i)
            t1,t2,t3=momentEstimate(each_quartet,upper)
            taus=(t1,t2,t3)
            push!(f,taus)
        end
        
        for each_tau in 1:length(f)
            if all(.>=(0), f[each_tau])
                positive=true
            else    
                positive=false
            end
        end

        if positive==true
            i=i+1
        elseif i==1
            upper=0.1
            break
        else
            upper=lower*((factor)^(i-1))
            break
        end
    end

    return lower,upper
end

"""
    golden_section_search(f,lower::Float64,upper::Float64,tolerance::Float64)

Using the loose interval for theta estimated from get_upper_lower_theta(), \\
here, we further tighten the interval for the feasible theta using the golden ratio.
"""
function golden_section_search(f,lower::Float64,upper::Float64,tolerance::Float64)
    goldenratio=2/(sqrt(5)+1)
    
    #Use the golden ratio to set the initial test points
    x1=upper - goldenratio*(upper-lower)
    x2=lower + goldenratio*(upper-lower)

    #evaluate function at the test points
    f1=f(x1)
    f2=f(x2)

    iteration = 0

    while (abs(upper-lower) > tolerance)
        iteration = iteration+1
        if f2 > f1
            upper=x2
            x2=x1
            f2=f1

            x1=upper-goldenratio*(upper-lower)
            f1=f(x1)
        else
            lower=x1
            x1=x2
            f1=f2

            x2=lower+goldenratio*(upper-lower)
            f2=f(x2)
        end
    end

    return lower,upper
end
