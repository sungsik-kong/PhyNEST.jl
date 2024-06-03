#written by Sungsik Kong 2021-2022
global const CoefSym1=Vector{Int}([1,-2,3,-2,6,-2,-8,-2,6])
global const CoefSym2=Vector{Int}([1,6,3,6,-2,-2,-8,-2,-2])
global const CoefSym3=Vector{Int}([1,2,-1,-2,2,2,-2,-2])
global const CoefASym1=Vector{Int}([3,-7,8,-8,10,-5,-10,-6,9,-12,18])
global const CoefASym2=Vector{Int}([3,5,-4,-8,-2,7,-10,-6,9,12,-6])
global const CoefASym3=Vector{Int}([3,10,1,2,5,2,4,-6,-3,-12,-6])
global const momestasymorder1=Vector{Int}([1,6,3,9,12,2,7,4,8,5,15,10,14,13,11])
global const momestasymorder2=Vector{Int}([1,6,3,9,12,10,4,7,14,13,15,2,8,5,11])
global const momestasymorder3=Vector{Int}([1,2,3,4,5,6,7,9,8,12,15,10,11,13,14])
global const momestasymorder4=Vector{Int}([1,10,6,4,14,3,7,9,13,12,15,2,11,8,5])
"""
methodofmomentestimator(each_quartet::quartets,theta::Float64)

For each quartet object, computes three branch lengths using the method of moment estimator
"""
function methodofmomentestimator(each_quartet::quartets,theta::Float64)
    quartet_is_symmetric=iszero(each_quartet.symtype)
    if quartet_is_symmetric momest=mom_est_sym(each_quartet.mspcountsNET,theta) 
        return momest
    elseif !(quartet_is_symmetric) momest=mom_est_asym(each_quartet.symtype,each_quartet.mspcountsNET,theta) 
        return momest
    end     
end

"""
    mom_est_sym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)

Given the site pattern frequencies, theta, and mu, this function computes three branch lengths using the method of moment estimator for a symmetric quartet.

## Mandatory arguments
- `site_pattern_frequnecies`      Vector of quartet site pattern frequencies
- `theta`     Effective population size parameter
- `alpha`     set dafault=4/3
"""
function mom_est_sym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)
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

    spprobS1=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyzxx]
    spprobS2=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,      pxyzw,pyzxx]
    
    myx1=4*(1+mu*2*theta)*sum(CoefSym1.*spprobS1)
    myx2=4*(1+mu*2*theta)*sum(CoefSym2.*spprobS1)
    myx3=4*(1+mu*2*theta)*sum(CoefSym3.*spprobS2)

    mytau1=-log(myx1^((1.0)/(2*mu)))/theta
    mytau2=-log(myx2^((1.0)/(2*mu)))/theta
    mytau3=-log(myx3^((1.0)/(2*mu)))/theta

    return (mytau1,mytau2,mytau3)
end

"""
    mom_est_asym(site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)

Given the site pattern frequencies, theta, and mu, computes three branch lengths using the method of moment estimator for the four types of asymmetric quartet.

## Mandatory arguments
- `type`      Type of the asymmetric quartet in interest. See above.
- `site_pattern_frequnecies`      Vector of quartet site pattern frequencies
- `theta`     Effective population size parameter
- `alpha`     set dafault=4/3
"""
function mom_est_asym(type::Integer,site_pattern_frequnecies::Array,theta::Float64; mu=(4/3)::Float64)

    function spprob(t)
        N=sum(site_pattern_frequnecies)
        pxxxx=site_pattern_frequnecies[t[1]]/(4*N)
        pxxxy=(site_pattern_frequnecies[t[2]]+site_pattern_frequnecies[t[3]])/(24*N)
        pxxyy=site_pattern_frequnecies[t[4]]/(12*N)
        pxxyz=site_pattern_frequnecies[t[5]]/(24*N)
        pxyxx=site_pattern_frequnecies[t[6]]/(12*N)
        pxyxy=(site_pattern_frequnecies[t[7]]+site_pattern_frequnecies[t[8]])/(24*N)
        pxyxz=(site_pattern_frequnecies[t[9]]+site_pattern_frequnecies[t[10]])/(48*N)
        pxyzw=site_pattern_frequnecies[t[11]]/(24*N)
        pyxxx=site_pattern_frequnecies[t[12]]/(12*N)
        pyxxz=(site_pattern_frequnecies[t[13]]+site_pattern_frequnecies[t[14]])/(48*N)
        pyzxx=site_pattern_frequnecies[t[15]]/(24*N)

        p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
        return p        
    end

    if type==1 spprobA=spprob(momestasymorder1)
    elseif type==2 spprobA=spprob(momestasymorder2)
    elseif type==3 spprobA=spprob(momestasymorder3)
    elseif type==4 spprobA=spprob(momestasymorder4)
    else error("There is no type $type.")
    end

    myx1=4/3*(1+mu*2*theta)*sum(CoefASym1.*spprobA)
    myx2=4/3*(1+mu*2*theta)*sum(CoefASym2.*spprobA)
    myx3=4/3*(1+mu*2*theta)*sum(CoefASym3.*spprobA)

    mytau1=-log(myx1^((1.0)/(2*mu)))/theta
    mytau2=-log(myx2^((1.0)/(2*mu)))/theta
    mytau3=-log(myx3^((1.0)/(2*mu)))/theta

    return (mytau1,mytau2,mytau3)
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


