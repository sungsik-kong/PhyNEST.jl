#written by Sungsik Kong 2021-2022
#Theta can be estimated, but if not provided, we use user provided theta value here (default=0.001)
"""
    momentEstimate(each_quartet::quartets,theta=0.001::Float64)

For each quartet object, computes three branch lengths using the method of moment estimator
"""
function momentEstimate(each_quartet::quartets,theta=0.001::Float64)
    sym_type=each_quartet.symtype
    if sym_type==0 
        t1,t2,t3=mom_est_sym(each_quartet.mspcountsNET,theta)
    else 
        t1,t2,t3=mom_est_asym(sym_type,each_quartet.mspcountsNET,theta) 
    end    
    return (t1,t2,t3)
end

"""
    mom_est_sym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)

Given the site pattern frequencies, theta, and mu, computes three branch lengths using \\
the method of moment estimator for a symmetric quartet.
"""
function mom_est_sym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)
    #mu=4/3
    CoefSym1=[1,-2,3,-2,6,-2,-8,-2,6]
    CoefSym2=[1,6,3,6,-2,-2,-8,-2,-2]
    CoefSym3=[1,2,-1,-2,2,2,-2,-2]

    spprobS1,spprobS2=mom_est_sp_prob_sym(site_pattern_frequnecies)
    
    myx1=4*(1+mu*2*theta)*sum(CoefSym1.*spprobS1)
    myx2=4*(1+mu*2*theta)*sum(CoefSym2.*spprobS1)
    myx3=4*(1+mu*2*theta)*sum(CoefSym3.*spprobS2)

    mytau1=-log(myx1^((1.0)/(2*mu)))/theta
    mytau2=-log(myx2^((1.0)/(2*mu)))/theta
    mytau3=-log(myx3^((1.0)/(2*mu)))/theta

    return mytau1,mytau2,mytau3
end

"""
    mom_est_sp_prob_sym(site_pattern_frequnecies::Array)

Computes site pattern probability for the observed site pattern frequencies for \\
computing branch lengths using the method of moment estimator for a symmetric quartet.
"""
function mom_est_sp_prob_sym(site_pattern_frequnecies::Array)
    N=sum(site_pattern_frequnecies)

    pxxxx=(site_pattern_frequnecies[1])/(4*N)
    pxxxy=(site_pattern_frequnecies[2]+site_pattern_frequnecies[3])/(24*N)
    pxxyy=(site_pattern_frequnecies[4])/(12*N)
    pxxyz=(site_pattern_frequnecies[5])/(24*N)
    pxyxx=(site_pattern_frequnecies[6]+site_pattern_frequnecies[10])/(24*N)
    pxyxy=(site_pattern_frequnecies[7]+site_pattern_frequnecies[9])/(24*N)
    pxyxz=(site_pattern_frequnecies[8]+site_pattern_frequnecies[11]+site_pattern_frequnecies[12]+site_pattern_frequnecies[13])/(96*N)
    pxyzw=(site_pattern_frequnecies[15])/(24*N)
    pyzxx=(site_pattern_frequnecies[14])/(24*N)

    p1=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyzxx]
    p2=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,      pxyzw,pyzxx]

    return p1,p2
end

"""
    mom_est_asym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)

Given the site pattern frequencies, theta, and mu, computes three branch lengths using \\
the method of moment estimator for the four types of asymmetric quartet.
"""
function mom_est_asym(type::Integer,site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)
    #mu=4/3
    CoefASym1=[3,-7,8,-8,10,-5,-10,-6,9,-12,18]
    CoefASym2=[3,5,-4,-8,-2,7,-10,-6,9,12,-6]
    CoefASym3=[3,10,1,2,5,2,4,-6,-3,-12,-6]

    if type==1
        spprobA=mom_est_asym_sp_prob_type1(site_pattern_frequnecies)
    elseif type==2
        spprobA=mom_est_asym_sp_prob_type2(site_pattern_frequnecies)
    elseif type==3
        spprobA=mom_est_asym_sp_prob_type3(site_pattern_frequnecies)
    elseif type==4
        spprobA=mom_est_asym_sp_prob_type4(site_pattern_frequnecies)
    else
        error("There is no type $type.")
    end

    myx1=4/3*(1+mu*2*theta)*sum(CoefASym1.*spprobA)
    myx2=4/3*(1+mu*2*theta)*sum(CoefASym2.*spprobA)
    myx3=4/3*(1+mu*2*theta)*sum(CoefASym3.*spprobA)

    mytau1=-log(myx1^((1.0)/(2*mu)))/theta
    mytau2=-log(myx2^((1.0)/(2*mu)))/theta
    mytau3=-log(myx3^((1.0)/(2*mu)))/theta

    return mytau1,mytau2,mytau3
end

"""
    mom_est_asym_sp_prob_type1(site_pattern_frequnecies::Array)

Computes site pattern probability for the observed site pattern frequencies for computing \\
branch lengths using the method of moment estimator for the asymmetric quartet type 1.
"""
function mom_est_asym_sp_prob_type1(site_pattern_frequnecies::Array)
    N=sum(site_pattern_frequnecies)

    pxxxx=site_pattern_frequnecies[1]/(4*N)
    pxxxy=(site_pattern_frequnecies[6]+site_pattern_frequnecies[3])/(24*N)
    pxxyy=site_pattern_frequnecies[9]/(12*N)
    pxxyz=site_pattern_frequnecies[12]/(24*N)
    pxyxx=site_pattern_frequnecies[2]/(12*N)
    pxyxy=(site_pattern_frequnecies[7]+site_pattern_frequnecies[4])/(24*N)
    pxyxz=(site_pattern_frequnecies[8]+site_pattern_frequnecies[5])/(48*N)
    pxyzw=site_pattern_frequnecies[15]/(24*N)
    pyxxx=site_pattern_frequnecies[10]/(12*N)
    pyxxz=(site_pattern_frequnecies[14]+site_pattern_frequnecies[13])/(48*N)
    pyzxx=site_pattern_frequnecies[11]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

"""
    mom_est_asym_sp_prob_type2(site_pattern_frequnecies::Array)

Computes site pattern probability for the observed site pattern frequencies for computing \\
branch lengths using the method of moment estimator for the asymmetric quartet type 2.
"""
function mom_est_asym_sp_prob_type2(site_pattern_frequnecies::Array)
    N=sum(site_pattern_frequnecies)

    pxxxx=site_pattern_frequnecies[1]/(4*N)
    pxxxy=(site_pattern_frequnecies[6]+site_pattern_frequnecies[3])/(24*N)
    pxxyy=site_pattern_frequnecies[9]/(12*N)
    pxxyz=site_pattern_frequnecies[12]/(24*N)
    pxyxx=site_pattern_frequnecies[10]/(12*N)
    pxyxy=(site_pattern_frequnecies[4]+site_pattern_frequnecies[7])/(24*N)
    pxyxz=(site_pattern_frequnecies[14]+site_pattern_frequnecies[13])/(48*N)
    pxyzw=site_pattern_frequnecies[15]/(24*N)
    pyxxx=site_pattern_frequnecies[2]/(12*N)
    pyxxz=(site_pattern_frequnecies[8]+site_pattern_frequnecies[5])/(48*N)
    pyzxx=site_pattern_frequnecies[11]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

"""
    mom_est_asym_sp_prob_type3(site_pattern_frequnecies::Array)

Computes site pattern probability for the observed site pattern frequencies for computing \\
branch lengths using the method of moment estimator for the asymmetric quartet type 3.
"""
function mom_est_asym_sp_prob_type3(site_pattern_frequnecies::Array)
    N=sum(site_pattern_frequnecies)

    pxxxx=site_pattern_frequnecies[1]/(4*N)
    pxxxy=(site_pattern_frequnecies[2]+site_pattern_frequnecies[3])/(24*N)
    pxxyy=site_pattern_frequnecies[4]/(12*N)
    pxxyz=site_pattern_frequnecies[5]/(24*N)
    pxyxx=site_pattern_frequnecies[6]/(12*N)
    pxyxy=(site_pattern_frequnecies[7]+site_pattern_frequnecies[9])/(24*N)
    pxyxz=(site_pattern_frequnecies[8]+site_pattern_frequnecies[12])/(48*N)
    pxyzw=site_pattern_frequnecies[15]/(24*N)
    pyxxx=site_pattern_frequnecies[10]/(12*N)
    pyxxz=(site_pattern_frequnecies[11]+site_pattern_frequnecies[13])/(48*N)
    pyzxx=site_pattern_frequnecies[14]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

"""
    mom_est_asym_sp_prob_type4(site_pattern_frequnecies::Array)

Computes site pattern probability for the observed site pattern frequencies for computing \\
branch lengths using the method of moment estimator for the asymmetric quartet type 4.
"""
function mom_est_asym_sp_prob_type4(site_pattern_frequnecies::Array)
    N=sum(site_pattern_frequnecies)

    pxxxx=site_pattern_frequnecies[1]/(4*N)
    pxxxy=(site_pattern_frequnecies[10]+site_pattern_frequnecies[6])/(24*N)
    pxxyy=site_pattern_frequnecies[4]/(12*N)
    pxxyz=site_pattern_frequnecies[14]/(24*N)
    pxyxx=site_pattern_frequnecies[3]/(12*N)
    pxyxy=(site_pattern_frequnecies[7]+site_pattern_frequnecies[9])/(24*N)
    pxyxz=(site_pattern_frequnecies[13]+site_pattern_frequnecies[12])/(48*N)
    pxyzw=site_pattern_frequnecies[15]/(24*N)
    pyxxx=site_pattern_frequnecies[2]/(12*N)
    pyxxz=(site_pattern_frequnecies[11]+site_pattern_frequnecies[8])/(48*N)
    pyzxx=site_pattern_frequnecies[5]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

"""
    get_average_moment_branch_length(N::Network)

Averages all estimated branch lengths using the method of moment estimator for each tau.
"""
function get_average_moment_branch_length(N::Network)
    all_quartets=N.quartet
    largest_ntau=0
    for each_quartet in all_quartets
        if each_quartet.ntau[1] > largest_ntau
            largest_ntau=each_quartet.ntau[1]
        elseif each_quartet.ntau[2] > largest_ntau
            largest_ntau=each_quartet.ntau[2]
        elseif each_quartet.ntau[3] > largest_ntau
            largest_ntau=each_quartet.ntau[3]
        end
    end

    average_moment_branch_length=Float64[]
    current_ntau=1
    while current_ntau<=largest_ntau
        current_moment_branch_length=Float64[]
        for each_quartet in all_quartets
            if each_quartet.ntau[1]==current_ntau
                push!(current_moment_branch_length,each_quartet.momestlength[1])
            elseif each_quartet.ntau[2]==current_ntau
                push!(current_moment_branch_length,each_quartet.momestlength[2])
            elseif each_quartet.ntau[3]==current_ntau
                push!(current_moment_branch_length,each_quartet.momestlength[3])
            end
        end
        push!(average_moment_branch_length,get_average(current_moment_branch_length))
        current_ntau+=1
    end
    return average_moment_branch_length
end

