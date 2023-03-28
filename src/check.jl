#Written by Sungsik Kong 2021-2022
#Last updated by Sungsik Kong 2023
function check_probabilities(t1::Float64,t2::Float64,t3::Float64,theta::Float64;alpha=(4/3))
    symmetric_probabilities=GetTrueProbsSymm(t1,t2,t3,theta,alpha)
    asymmetric_probabilities_type1=GetTrueProbsAsymmTypes(1,t1,t2,t3,theta,alpha)
    asymmetric_probabilities_type2=GetTrueProbsAsymmTypes(2,t1,t2,t3,theta,alpha)
    asymmetric_probabilities_type3=GetTrueProbsAsymmTypes(3,t1,t2,t3,theta,alpha)
    asymmetric_probabilities_type4=GetTrueProbsAsymmTypes(4,t1,t2,t3,theta,alpha)

    weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]

    symmetric_probabilities=(weights .* symmetric_probabilities)
    asymmetric_probabilities_type1=(weights .* asymmetric_probabilities_type1)
    asymmetric_probabilities_type2=(weights .* asymmetric_probabilities_type2)
    asymmetric_probabilities_type3=(weights .* asymmetric_probabilities_type3)
    asymmetric_probabilities_type4=(weights .* asymmetric_probabilities_type4)

    println(sum(symmetric_probabilities))
    println(sum(asymmetric_probabilities_type1))
    println(sum(asymmetric_probabilities_type2))
    println(sum(asymmetric_probabilities_type3))
    println(sum(asymmetric_probabilities_type4))
end

function check_mom_estimates_and_theta(t1::Float64,t2::Float64,t3::Float64,theta::Float64;alpha=(4/3)::Float64,n=10000000000::Int64)
    tree_type_0=readTopology("((i,j),(k,l));")
    tree_type_1=readTopology("(i,((j,k),l));")
    tree_type_2=readTopology("((i,(j,k)),l);") 
    tree_type_3=readTopology("(i,(j,(k,l)));")
    tree_type_4=readTopology("(((i,j),k),l);")

    N0=get_quartets(tree_type_0)
    N1=get_quartets(tree_type_1)
    N2=get_quartets(tree_type_2)
    N3=get_quartets(tree_type_3)
    N4=get_quartets(tree_type_4)

    N0_q=N0.quartet[1]
    N1_q=N1.quartet[1]
    N2_q=N2.quartet[1]
    N3_q=N3.quartet[1]
    N4_q=N4.quartet[1]

    N0_q.mspcountsNET=sim_sp_counts(0,t1,t2,t3,theta,alpha,n)
    N1_q.mspcountsNET=sim_sp_counts(1,t1,t2,t3,theta,alpha,n)
    N2_q.mspcountsNET=sim_sp_counts(2,t1,t2,t3,theta,alpha,n)
    N3_q.mspcountsNET=sim_sp_counts(3,t1,t2,t3,theta,alpha,n)
    N4_q.mspcountsNET=sim_sp_counts(4,t1,t2,t3,theta,alpha,n)

    mom_N0_t1,mom_N0_t2,mom_N0_t3=momentEstimate(N0_q,theta)
    mom_N1_t1,mom_N1_t2,mom_N1_t3=momentEstimate(N1_q,theta)
    mom_N2_t1,mom_N2_t2,mom_N2_t3=momentEstimate(N2_q,theta)
    mom_N3_t1,mom_N3_t2,mom_N3_t3=momentEstimate(N3_q,theta)
    mom_N4_t1,mom_N4_t2,mom_N4_t3=momentEstimate(N4_q,theta)

    theta_N0=get_start_theta(N0)
    theta_N1=get_start_theta(N1)
    theta_N2=get_start_theta(N2)
    theta_N3=get_start_theta(N3)
    theta_N4=get_start_theta(N4)

    println([N0_q.symtype,mom_N0_t1,mom_N0_t2,mom_N0_t3,theta_N0])
    println([N1_q.symtype,mom_N1_t1,mom_N1_t2,mom_N1_t3,theta_N1])
    println([N2_q.symtype,mom_N2_t1,mom_N2_t2,mom_N2_t3,theta_N2])
    println([N3_q.symtype,mom_N3_t1,mom_N3_t2,mom_N3_t3,theta_N3])
    println([N4_q.symtype,mom_N4_t1,mom_N4_t2,mom_N4_t3,theta_N4])
end