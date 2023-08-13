#Written by Sungsik Kong 2021-2022
#Last updated by Sungsik Kong 2023

"""
    GetTrueProbsSymm(myt1::Float64,
                    myt2::Float64,
                    myt3::Float64,
                    theta::Float64)
    GetTrueProbsSymm(myt1::Float64,
                    myt2::Float64,
                    myt3::Float64,
                    theta::Float64,
                    alpha::Float64)

Computes true site pattern probabilities for the symmetric quartet tree: `((1,2),(3,4));`. 
Three speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified. 
See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.

## Input
-`myt1`      Speciation time for the common ancestor of species 1 and 2 in coalescent unit
-`myt2`      Speciation time for the common ancestor of species 3 and 4 in coalescent unit
-`myt3`      Root node age in coalescent unit
-`theta`     Effective population size parameter
-`alpha`     (dafault=4/3)
"""
GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)=GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,4/3)
function GetTrueProbsSymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)

    t1=myt1*theta
    t2=myt2*theta
    t3=myt3*theta
    t=2*theta
    m=alpha
    
    #compute the C matrix
    cmat=zeros(Float64,9,9)
    
    #row 1
    for i in 1:9 cmat[1,i] = 1/256 end
    #row 2
    cmat[2,1]=cmat[2,3]=cmat[2,5]=cmat[2,7] = 3/(256*(1+m*t))
    cmat[2,2]=cmat[2,4]=cmat[2,6]=cmat[2,8]=cmat[2,9]=-1/(256*(1+m*t))
    # row 3
    cmat[3,1]=cmat[3,2]=cmat[3,5]=cmat[3,6]=3/(256*(1+m*t))
    cmat[3,3]=cmat[3,4]=cmat[3,7]=cmat[3,8]=cmat[3,9]=-1/(256*(1+m*t))
    # row 4
    cmat[4,1]=cmat[4,5] = 9/(256*(1+m*t)^2)
    cmat[4,2]=cmat[4,3]=cmat[4,6]=cmat[4,7] = -3/(256*(1+m*t)^2)
    cmat[4,4]=cmat[4,8]=cmat[4,9] = 1/(256*(1+m*t)^2)
    # row 5
    cmat[5,1]=12/(256*(1+m*t))
    cmat[5,2]=cmat[5,3]=cmat[5,4]=4/(256*(1+m*t))
    cmat[5,5]=cmat[5,6]=cmat[5,7]=cmat[5,9]=-4/(256*(1+m*t))
    cmat[5,8]=0
    # row 6
    cmat[6,1]=24/(256*(1+m*t)*(2+m*t))
    cmat[6,2]=cmat[6,4]=cmat[6,5]=cmat[6,7]=-8/(256*(1+m*t)*(2+m*t))
    cmat[6,3]=8/(256*(1+m*t)*(2+m*t))
    cmat[6,6]=cmat[6,9]=8/(256*(1+m*t)*(2+m*t))
    cmat[6,8]=0
    # row 7
    cmat[7,1]=24/(256*(1+m*t)*(2+m*t))
    cmat[7,2]=8/(256*(1+m*t)*(2+m*t))
    cmat[7,3]=cmat[7,4]=cmat[7,5]=cmat[7,6]=-8/(256*(1+m*t)*(2+m*t))
    cmat[7,7]=cmat[7,9]=8/(256*(1+m*t)*(2+m*t))
    cmat[7,8]=0
    # row 8
    cmat[8,1]=48/(256*(1+m*t)*(2+m*t)^2)
    cmat[8,2]=cmat[8,3]=cmat[8,5]=-16/(256*(1+m*t)*(2+m*t)^2)
    cmat[8,4]=cmat[8,6]=cmat[8,7]=16/(256*(1+m*t)*(2+m*t)^2)
    cmat[8,8]=0
    cmat[8,9]=-16/(256*(1+m*t)*(2+m*t)^2)
    # row 9
    cmat[9,1]=6*m*t*(4+m*t)*(4+3*m*t)/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,2]=cmat[9,3]=-2*m*t*(4+m*t)*(4+3*m*t)/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,4]=m*t*(32+40*m*t+10*(m^2)*(t^2))/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,5]=2*m*t*(4+m*t)^2/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,6]=cmat[9,7]=2*(m^2)*(t^2)*(4+m*t)/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,8]=-(m^2)*(t^2)*(4+2*m*t)/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    cmat[9,9]=2*(m^3)*(t^3)/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
        
    beta=zeros(Float64,9,1)
    beta[1,1]=1
    beta[2,1]=exp(-2*m*t1)
    beta[3,1]=exp(-2*m*t2)
    beta[4,1]=exp(-2*m*t1)*exp(-2*m*t2)
    beta[5,1]=exp(-2*m*t3)
    beta[6,1]=exp(-m*t1)*exp(-2*m*t3)
    beta[7,1]=exp(-m*t2)*exp(-2*m*t3)
    beta[8,1]=exp(-m*t1)*exp(-m*t2)*exp(-2*m*t3)
    beta[9,1]=exp(2*t1/t)*exp(2*t2/t)*exp(-4*t3*(m+1/t))	    

    p=(transpose(cmat))*beta

    p1=zeros(Float64,15)
    p1[1]=p[1,1] ##xxxx
    p1[2]=p[2,1] ##xxxy
    p1[3]=p[2,1] ##xxyx
    p1[4]=p[5,1] ##xxyy
    p1[5]=p[6,1] ##xxyz
    p1[6]=p[3,1] ##xyxx
    p1[7]=p[4,1] ##xyxy
    p1[8]=p[8,1] ##xyxz
    p1[9]=p[4,1] ##xyyx
    p1[10]=p[3,1] ##yxxx
    p1[11]=p[8,1] ##xyyz
    p1[12]=p[8,1] ##zxyz
    p1[13]=p[8,1] ##yxzx
    p1[14]=p[7,1] ##yzxx
    p1[15]=p[9,1] ##xyzw
    
    return p1
end


"""
    GetTrueProbsAsymm(myt1::Float64,
                    myt2::Float64,
                    myt3::Float64,
                    theta::Float64)
    GetTrueProbsAsymm(myt1::Float64,
                    myt2::Float64,
                    myt3::Float64,
                    theta::Float64,
                    alpha::Float64)

Computes true site pattern probabilities for the asymmetric quartet tree: (1,(2,(3,4)));.
Three speciation times (node ages) in coalescent unit and theta must be provided; alpha is assumed to be 4/3 if unspecified. 
See the manuscript and/or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.

## Input
-`myt1`      Speciation time for the common ancestor of species 3 and 4 in coalescent unit
-`myt2`      Speciation time for the common ancestor of species 2 and (3,4) in coalescent unit
-`myt3`      Root node age in coalescent unit
-`theta`     Effective population size parameter
-`alpha`     (dafault=4/3)
"""
GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)=GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,4/3)
function GetTrueProbsAsymm(myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)

    t1 = myt1*theta
    t2 = myt2*theta  
    t3 = myt3*theta
    t = 2*theta
    m = alpha #mu for JC69
    
    #compute the C matrix
    cmat = zeros(Float64,11,10)
    
    # row 1
    cmat[1,1] = 1/256
    cmat[1,2] = 3/((256)*(1+m*t))
    cmat[1,3] = 6/(256*(1+m*t))
    cmat[1,4] = 12/(256*(1+m*t)*(2+m*t))
    cmat[1,5] = 9/(256*(1+m*t))
    cmat[1,6] = 12/(256*(1+m*t)*(2+m*t))#6/(256*(1+m*t)*(2+m*t))#12/(256*(1+m*t)*(2+m*t))
    cmat[1,7] = 9/(256*(1+m*t)^2)
    cmat[1,8] = 24/(256*(1+m*t)*(2+m*t))
    cmat[1,9] = 48/(256*(1+m*t)*(2+m*t)^2)
    cmat[1,10] =(6*m*t*(4+m*t)*(4+3*m*t))/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 2
    cmat[2,1] = 1/256
    cmat[2,2] = -1/(256*(1+m*t))
    cmat[2,3] = 2/(256*(1+m*t))
    cmat[2,4] = -4/(256*(1+m*t)*(2+m*t))
    cmat[2,5] = 5/(256*(1+m*t))
    cmat[2,6] = -4/(256*(1+m*t)*(2+m*t))
    cmat[2,7] = -3/(256*(1+m*t)^2)
    cmat[2,8] = 8/(256*(1+m*t)*(2+m*t))
    cmat[2,9] = -16/(256*(1+m*t)*(2+m*t)^2)
    cmat[2,10] = -(2*m*t*(4+m*t)*(4+3*m*t))/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    
    # row 3
    cmat[3,1] = 1/256
    cmat[3,2] = 3/(256*(1+m*t))
    cmat[3,3] = -2/(256*(1+m*t))
    cmat[3,4] = -4/(256*(1+m*t)*(2+m*t))
    cmat[3,5] = 5/(256*(1+m*t))
    cmat[3,6] = 12/(256*(1+m*t)*(2+m*t))
    cmat[3,7] = -3/(256*(1+m*t)^2)
    cmat[3,8] = -8/(256*(1+m*t)*(2+m*t))
    cmat[3,9] = -16/(256*(1+m*t)*(2+m*t)^2)
    cmat[3,10] = -(2*m*t*(4+m*t)*(4+3*m*t))/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))

    # row 4
    cmat[4,1] = 1/256
    cmat[4,2] = 3/(256*(1+m*t))
    cmat[4,3] = 6/(256*(1+m*t))
    cmat[4,4] = 12/(256*(1+m*t)*(2+m*t))
    cmat[4,5] = -3/(256*(1+m*t))
    cmat[4,6] = -4/(256*(1+m*t)*(2+m*t))
    cmat[4,7] = -3/(256*(1+m*t)^2)
    cmat[4,8] = -8/(256*(1+m*t)*(2+m*t))
    cmat[4,9] = -16/(256*(1+m*t)*(2+m*t)^2)
    cmat[4,10] = -(2*m*t*(4+m*t)*(4+3*m*t))/(256*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    #row 5
    cmat[5,1] = 1/256
    cmat[5,2] = 3/((256)*(1+m*t))
    cmat[5,3] = -2/((256)*(1+m*t))
    cmat[5,4] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[5,5] = 1/((256)*(1+m*t))
    cmat[5,6] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[5,7] = 9/((256)*(1+m*t)^2)
    cmat[5,8] = -8/((256)*(1+m*t)*(2+m*t))
    cmat[5,9] = -16/((256)*(1+m*t)*(2+m*t)^2)
    cmat[5,10] = (2*m*t*(4+m*t)^2)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 6
    cmat[6,1] = 1/256
    cmat[6,2] = -1/((256)*(1+m*t))
    cmat[6,3] = 2/((256)*(1+m*t))
    cmat[6,4] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[6,5] = 1/((256)*(1+m*t))
    cmat[6,6] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[6,7] = 1/((256)*(1+m*t)^2)
    cmat[6,8] = -8/((256)*(1+m*t)*(2+m*t))
    cmat[6,9] = 16/((256)*(1+m*t)*(2+m*t)^2)
    cmat[6,10] = m*t*(2*16+40*m*t+(10)*(m^2)*(t^2))/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 7
    cmat[7,1] = 1/256
    cmat[7,2] = -1/((256)*(1+m*t))
    cmat[7,3] = -2/((256)*(1+m*t))
    cmat[7,4] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[7,5] = 1/((256)*(1+m*t))
    cmat[7,6] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[7,7] = -3/((256)*(1+m*t)^2)#-1/((256)*(1+m*t)^2)#-3/((256)*(1+m*t)^2)
    cmat[7,8] = -8/((256)*(1+m*t)*(2+m*t))
    cmat[7,9] = 16/((256)*(1+m*t)*(2+m*t)^2)
    cmat[7,10] = 2*(m^2)*(t^2)*(4+m*t)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 8
    cmat[8,1] = 1/256
    cmat[8,2] = 3/((256)*(1+m*t))
    cmat[8,3] = -2/((256)*(1+m*t))
    cmat[8,4] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[8,5] = -3/((256)*(1+m*t))
    cmat[8,6] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[8,7] = -3/((256)*(1+m*t)^2)
    cmat[8,8] = 8/((256)*(1+m*t)*(2+m*t))
    cmat[8,9] = 16/((256)*(1+m*t)*(2+m*t)^2)#32/((256)*(1+m*t)*(2+m*t)^2)#16/((256)*(1+m*t)*(2+m*t)^2)
    cmat[8,10] = 2*(m^2)*(t^2)*(4+m*t)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 9
    cmat[9,1] = 1/256
    cmat[9,2] = -1/((256)*(1+m*t))
    cmat[9,3] = -2/((256)*(1+m*t))
    cmat[9,4] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[9,5] = 1/((256)*(1+m*t))
    cmat[9,6] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[9,7] = 1/((256)*(1+m*t)^2)
    cmat[9,8] = 0
    cmat[9,9] = 0
    cmat[9,10] = -(m^2)*(t^2)*(4+2*m*t)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 10
    cmat[10,1] = 1/256
    cmat[10,2] = -1/((256)*(1+m*t))
    cmat[10,3] = 2/((256)*(1+m*t))
    cmat[10,4] = -4/((256)*(1+m*t)*(2+m*t))
    cmat[10,5] = -3/((256)*(1+m*t))
    cmat[10,6] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[10,7] = 1/((256)*(1+m*t)^2)
    cmat[10,8] = 0
    cmat[10,9] = 0
    cmat[10,10] = -(m^2)*(t^2)*(4+2*m*t)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))
    # row 11
    cmat[11,1] = 1/256
    cmat[11,2] = -1/((256)*(1+m*t))
    cmat[11,3] = -2/((256)*(1+m*t))
    cmat[11,4] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[11,5] = -3/((256)*(1+m*t))
    cmat[11,6] = 4/((256)*(1+m*t)*(2+m*t))
    cmat[11,7] = 1/((256)*(1+m*t)^2)
    cmat[11,8] = 8/((256)*(1+m*t)*(2+m*t))
    cmat[11,9] = -16/((256)*(1+m*t)*(2+m*t)^2)
    cmat[11,10] = 2*(m^3)*(t^3)/((256)*((1+m*t)^2)*((2+m*t)^2)*(3+m*t))

    #get the beta vector
    beta = zeros(Float64,10,1)
    beta[1,1]=1
    beta[2,1]=exp(-2*m*t1)
    beta[3,1]=exp(-2*m*t2)
    beta[4,1]=exp(-m*t1)*exp(-2*m*t2)
    beta[5,1]=exp(-2*m*t3)
    beta[6,1]=exp(-m*t1)*exp(-2*m*t3)
    beta[7,1]=exp(-2*m*t1)*exp(-2*m*t3)
    beta[8,1]=exp(-m*t2)*exp(-2*m*t3)
    beta[9,1]=exp(-m*t1)*exp(-m*t2)*exp(-2*m*t3)
    beta[10,1]=exp((2/t)*(t1-t2))*exp(-2*m*(t2+t3))

    p=cmat*beta

    p1=zeros(Float64,15)
    p1[1]=p[1,1] ##xxxx
    p1[2]=p[2,1] ##xxxy
    p1[3]=p[2,1] ##xxyx
    p1[4]=p[5,1] ##xxyy
    p1[5]=p[7,1] ##xxyz
    p1[6]=p[3,1] ##xyxx
    p1[7]=p[6,1] ##xyxy
    p1[8]=p[9,1] ##xyxz
    p1[9]=p[6,1] ##xyyx
    p1[10]=p[4,1] ##yxxx
    p1[11]=p[10,1] ##xyyz
    p1[12]=p[9,1] ##zxyz
    p1[13]=p[10,1] ##yxzx
    p1[14]=p[8,1] ##yzxx
    p1[15]=p[11,1] ##xyzw
   
    return p1

end











































"""
    GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)
    GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)

Computes true site pattern probabilities for any of the four the asymmetric quartet trees: \\

- Type 1: (i,((j,k),l));
- Type 2: ((i,(j,k)),l); 
- Type 3: (i,(j,(k,l))); 
- Type 4: (((i,j),k),l).\\
Three speciation times (or node ages) in coalescent unit and theta must be provided. \\
Alhpa is not an essential argument for the function and if not provided, it isssumed to be 4/3 by default. \\
See manuscript or Chifman and Kubatko (2015)[10.1016/j.jtbi.2015.03.006] for more information.\\

## Input
`myt1`      Speciation time for the common ancestor of species 3 and 4\\
`myt2`      Speciation time for the common ancestor of species 2, 3 and 4\\
`myt3`      Root node age\\
`theta`     Effective population size parameter\\
`alpha`     4/3 by default
"""
GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64)=GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,4/3)
function GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
    p=GetTrueProbsAsymm(myt1,myt2,myt3,theta,alpha)
    newp=zeros(Float64,15)#1432
    if type==1
        i=[1,6,3,9,12,2,7,8,4,10,14,5,13,11,15]
    elseif type==2
        i=[1,10,2,9,11,3,7,13,4,6,14,5,8,12,15]
    elseif type==3
        i=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
    elseif type==4
        i=[1,10,6,4,14,3,7,13,9,2,11,12,8,5,15]
    else
        error("There is no asymmetric quartet type $type. It should be between 1 and 4.")
    end

    for n in 1:15
        newp[n]=p[i[n]]
    end
    return newp 
end

#shame to the bad function naming skill...
function GetTrueProbsNetTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64)
    GetTrueProbsAsymmTypes(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64) end

""" 
    sim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64)
    sim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,length::Integer)

This function simulates the 15 site pattern probabilities for five quartet tree topologies, \\
one symmetric and four asymmetric quartets: 

- Type 0: ((i,j),(k,l));
- Type 1: (i,((j,k),l));
- Type 2: ((i,(j,k)),l); 
- Type 3: (i,(j,(k,l))); 
- Type 4: (((i,j),k),l).\\

## Input
`type`   Specify the type of a quartet (use integer)\\
`myt1`   Speciation time for the most recent internal tree node. Common ancestor of 1 and 2 in the symmetric case\\
`myt2`   Speciation time for the internal tree node closer to the root\\
`myt3`   Root node age\\
`theta`  Effective population size parameter (default=0.01)\\
`alpha`  Alpha (default=4/3)\\
`length` Sequence lengths (default=1000000)
"""
sim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64)=sim_sp_counts(type,myt1,myt2,myt3,0.01,4/3,1000000)
function sim_sp_counts(type::Integer,myt1::Float64,myt2::Float64,myt3::Float64,theta::Float64,alpha::Float64,n::Integer)
    weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
    if type==0
        p=GetTrueProbsSymm(myt1,myt2,myt3,theta,alpha)
    else
        p=GetTrueProbsNetTypes(type,myt1,myt2,myt3,theta,alpha)
    end
    prob=p.*weights
    sim=rand(Multinomial(n,prob))
    return sim
end
