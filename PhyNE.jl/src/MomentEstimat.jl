#written by Sungsik Kong 2021-2022

function momentEstimat(q::nquartets,theta::Float64)
    type=q.symtype[1]
    
    if type==0 t1,t2,t3=momEstSym(q.mspcountsNET[1],theta)
    else t1,t2,t3=MOMEstAsym(type,q.mspcountsNET[1],theta) 
    end    
            
    return t1,t2,t3
end

function momentEstimat(type::Integer,spcounts::Array,theta::Float64)
    
    if type==0 t1,t2,t3=momEstSym(spcounts,theta)
    else t1,t2,t3=MOMEstAsym(type,spcounts,theta) 
    end    
            
    return t1,t2,t3
end

function momEstSym(q::nquartets,thetaa::Float64)
    s=q.mspcounts[1]

    t1,t2,t3=momEstSym(s,thetaa)

    return t1,t2,t3
end

function momEstSym(sitefreq::Array,theta::Float64)
    mu=4/3
    CoefSym1=Array{Int64,1}([1,-2,3,-2,6,-2,-8,-2,6])
    CoefSym2=Array{Int64,1}([1,6,3,6,-2,-2,-8,-2,-2])
    CoefSym3=Array{Int64,1}([1,2,-1,-2,2,2,-2,-2])

    spprobS1,spprobS2=spProbSym(sitefreq)
    
    myx1=4*(1+mu*2*theta)*sum(CoefSym1.*spprobS1)
    myx2=4*(1+mu*2*theta)*sum(CoefSym2.*spprobS1)
    myx3=4*(1+mu*2*theta)*sum(CoefSym3.*spprobS2)

    mytau1=-log(myx1^((1.0)/(2*mu)))/theta
    mytau2=-log(myx2^((1.0)/(2*mu)))/theta
    mytau3=-log(myx3^((1.0)/(2*mu)))/theta

    return mytau1,mytau2,mytau3
end

function spProbSym(mspcounts::Array)
    N=sum(mspcounts)

    pxxxx=(mspcounts[1])/(4*N)
    pxxxy=(mspcounts[2]+mspcounts[3])/(24*N)
    pxxyy=(mspcounts[4])/(12*N)
    pxxyz=(mspcounts[5])/(24*N)
    pxyxx=(mspcounts[6]+mspcounts[10])/(24*N)
    pxyxy=(mspcounts[7]+mspcounts[9])/(24*N)
    pxyxz=(mspcounts[8]+mspcounts[11]+mspcounts[12]+mspcounts[13])/(96*N)
    pxyzw=(mspcounts[15])/(24*N)
    pyzxx=(mspcounts[14])/(24*N)

    p1=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyzxx]
    p2=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,      pxyzw,pyzxx]

    return p1,p2
end

function MOMEstAsym(type::Integer,spcounts::Array,theta::Float64)

    mu=4/3
    CoefASym1=Array{Int64,1}([3,-7,8,-8,10,-5,-10,-6,9,-12,18])
    CoefASym2=Array{Int64,1}([3,5,-4,-8,-2,7,-10,-6,9,12,-6])
    CoefASym3=Array{Int64,1}([3,10,1,2,5,2,4,-6,-3,-12,-6])

    if type==1
        spprobA=MOMspProbAsym1(spcounts)
    elseif type==2
        spprobA=MOMspProbAsym2(spcounts)
    elseif type==3
        spprobA=MOMspProbAsym3(spcounts)
    elseif type==4
        spprobA=MOMspProbAsym4(spcounts)
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

function MOMspProbAsym1(mspcounts::Array)
    N=sum(mspcounts)

    pxxxx=mspcounts[1]/(4*N)
    pxxxy=(mspcounts[6]+mspcounts[3])/(24*N)
    pxxyy=mspcounts[9]/(12*N)
    pxxyz=mspcounts[12]/(24*N)
    pxyxx=mspcounts[2]/(12*N)
    pxyxy=(mspcounts[7]+mspcounts[4])/(24*N)
    pxyxz=(mspcounts[8]+mspcounts[5])/(48*N)
    pxyzw=mspcounts[15]/(24*N)
    pyxxx=mspcounts[10]/(12*N)
    pyxxz=(mspcounts[14]+mspcounts[13])/(48*N)
    pyzxx=mspcounts[11]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

function MOMspProbAsym2(mspcounts::Array)
    N=sum(mspcounts)

    pxxxx=mspcounts[1]/(4*N)
    pxxxy=(mspcounts[6]+mspcounts[3])/(24*N)
    pxxyy=mspcounts[9]/(12*N)
    pxxyz=mspcounts[12]/(24*N)
    pxyxx=mspcounts[10]/(12*N)
    pxyxy=(mspcounts[4]+mspcounts[7])/(24*N)
    pxyxz=(mspcounts[14]+mspcounts[13])/(48*N)
    pxyzw=mspcounts[15]/(24*N)
    pyxxx=mspcounts[2]/(12*N)
    pyxxz=(mspcounts[8]+mspcounts[5])/(48*N)
    pyzxx=mspcounts[11]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

function MOMspProbAsym3(mspcounts::Array)
    N=sum(mspcounts)

    pxxxx=mspcounts[1]/(4*N)
    pxxxy=(mspcounts[2]+mspcounts[3])/(24*N)
    pxxyy=mspcounts[4]/(12*N)
    pxxyz=mspcounts[5]/(24*N)
    pxyxx=mspcounts[6]/(12*N)
    pxyxy=(mspcounts[7]+mspcounts[9])/(24*N)
    pxyxz=(mspcounts[8]+mspcounts[12])/(48*N)
    pxyzw=mspcounts[15]/(24*N)
    pyxxx=mspcounts[10]/(12*N)
    pyxxz=(mspcounts[11]+mspcounts[13])/(48*N)
    pyzxx=mspcounts[14]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end

function MOMspProbAsym4(mspcounts::Array)
    N=sum(mspcounts)

    pxxxx=mspcounts[1]/(4*N)
    pxxxy=(mspcounts[10]+mspcounts[6])/(24*N)
    pxxyy=mspcounts[4]/(12*N)
    pxxyz=mspcounts[14]/(24*N)
    pxyxx=mspcounts[3]/(12*N)
    pxyxy=(mspcounts[7]+mspcounts[9])/(24*N)
    pxyxz=(mspcounts[13]+mspcounts[12])/(48*N)
    pxyzw=mspcounts[15]/(24*N)
    pyxxx=mspcounts[2]/(12*N)
    pyxxz=(mspcounts[11]+mspcounts[8])/(48*N)
    pyzxx=mspcounts[5]/(24*N)

    p=[pxxxx,pxxxy,pxxyy,pxxyz,pxyxx,pxyxy,pxyxz,pxyzw,pyxxx,pyxxz,pyzxx]
    
    return p
end