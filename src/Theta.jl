#written by Sungsik Kong 2021-2022

function upperLower(q::Array{Nquartets, 1};lbound=0.00001::Float64, factor=2.0::Float64)
    
    i=1
    lower=lbound
#   upper=lower*10^i
    upper=lower*((factor)^i)
    positive=true

    while (positive=true)
        f=[]
        for eachT in q
            for eachQ in eachT.nquartet
#               upper=lower*10^i
                upper=lower*((factor)^i)
                t1,t2,t3=momentEstimat(eachQ,upper)
                taus=[t1,t2,t3]
                push!(f,taus)
            end
        end
        
        for g in 1:length(f)
            if all(.>=(0), f[g])
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
#           upper=lower*10^(i-1)
            upper=lower*((factor)^(i-1))
            break
        end
    end

    return lower,upper
end

function upperLower(p::Phylip,net::HybridNetwork;lbound=0.00001::Float64,factor=2.0::Float64)
    q=extractNQuartets(net,p)
    lower,upper=upperLower(q,lbound=lbound,factor=factor)
    return lower,upper
end

function goldenSectionSearch(f,lower::Float64,upper::Float64,tolerance::Float64)
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

pushMomEstlength(newT::Array{Nquartets, 1},theta::Float64)=pushMomEstlength(newT,theta,false)
pushMomEstlength(net::HybridNetwork,p::Phylip,theta::Float64)=pushMomEstlength(net,p,theta,true)
function pushMomEstlength(newT::Array{Nquartets, 1},theta::Float64,check::Bool)
    for eachT in newT
        for eachQ in eachT.nquartet
            t1,t2,t3=momentEstimat(eachQ,theta)
            push!(eachQ.momestlength,[t1,t2,t3])
            if(check) println(eachQ.momestlength) end
        end
    end

    return newT
end

function pushMomEstlength(net::HybridNetwork,p::Phylip,theta::Float64,check::Bool)
    q=extractNQuartets(net,p)
    pushMomEstlength(q,theta,check)
end

function getmomBLAve(newT::Array{Nquartets, 1},net::HybridNetwork)
    numTaus=net.numNodes-net.numTaxa-net.numHybrids
    for n in 1:length(net.node)
        if net.node[n].number<0
            if tauNum(net.node[n].number) > numTaus
                numTaus=tauNum(net.node[n].number)
            else continue
            end
        else continue
        end
    end

    ave=Array[]
    for eachT in newT
        for eachQ in eachT.nquartet
            taus=zeros(Float64,(numTaus,1))
            i=eachQ.ntau[1][1]
            j=eachQ.ntau[1][2]
            k=eachQ.ntau[1][3]
            
            taus[i]=eachQ.momestlength[1][1]
            taus[j]=eachQ.momestlength[1][2]
            taus[k]=eachQ.momestlength[1][3]
            push!(ave,taus)               
        end
    end
    #then get a average in the end
    avef=Float64[]
    for eacht in 1:numTaus
        avee=Float64[]
        for eachQ in ave
            if !iszero(eachQ[eacht]); push!(avee,eachQ[eacht]) end
        end
        push!(avef,mean(avee))
    end
    return avef
end

function getmomBLAve(net::HybridNetwork,p::Phylip)
    q=extractNQuartets(net,p)
    q=pushMomEstlength(x,p,theta,false)
    avef=getmomBLAve(q,net)
    println(avef)
    return avef
end

#check
#ave=getmomBLAve(x,p)

pushMomBL(newT::Array{Nquartets, 1},avef::Array)=pushMomBL(newT,avef,false)
function pushMomBL(newT::Array{Nquartets, 1},avef::Array,check::Bool)
    for eachT in newT
        for eachQ in eachT.nquartet
            push!(eachQ.mombl,[avef[eachQ.ntau[1][1]],avef[eachQ.ntau[1][2]],avef[eachQ.ntau[1][3]]])
            if(check) println(eachQ.mombl) end
        end
    end

    return newT
end

pushMomBL(newT::Array{Nquartets, 1},net::HybridNetwork)=pushMomBL(newT,net,false)
function pushMomBL(newT::Array{Nquartets, 1},net::HybridNetwork,check::Bool)
    avef=getmomBLAve(newT,net)
    newT=pushMomBL(newT,avef,check)
    return newT
end

pushMomBL(net::HybridNetwork,p::Phylip)=pushMomBL(net,p,true)
function pushMomBL(net::HybridNetwork,p::Phylip,check::Bool)
    q=extractNQuartets(net,p)
    q=pushMomEstlength(x,p,theta,false)
    pushMomBL(q,net,check)
end

function GetTrueProbsTh(q::nquartets,theta::Float64,alpha::Float64,gamma::Float64)
    
    type=q.symtype[1]
    t1=q.mombl[1][1]
    t2=q.mombl[1][2]
    t3=q.mombl[1][3]
    
    if type==0 p=gamma*GetTrueProbsSymm(t1,t2,t3,theta,alpha)
            else p=gamma*GetTrueProbsNetTypes(type,t1,t2,t3,theta,alpha)
    end

    return p
    
end

function getLikelihoodTh(q::Array{Nquartets, 1},theta1::Float64,alpha::Float64)
#    printQuarts(q)
    weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
    for eachT in q
        gam=1/length(q)

        for eachQ in eachT.nquartet

            myprobvec1 = GetTrueProbsTh(eachQ,theta1,alpha,gam)

            probs = (weights .* myprobvec1)    
            sitefreq=eachQ.mspcountsNET[1]
            if (sum(ismissing.(probs))!=0) 
                eachQ.logLik=-10^308
            elseif (sum(any(t->t<=0,probs))==0)
                lik=BigFloat(1.0)
                for i in 1:15
                    lik*=BigFloat(probs[i])^sitefreq[i]
                end
                eachQ.logLik=BigFloat(lik)
            else 
                eachQ.logLik=-10^308
            end 
        end
    end
    return q
end

function getLogPseudoLik(newT::Array{Nquartets, 1})
    Pseudolikelihood=BigFloat(0)
    for eachT in newT
        for eachQ in eachT.nquartet
            Pseudolikelihood+=eachQ.logLik
        end
    end

    y=-1.0*(log(Pseudolikelihood))

    return Float64(y)
end

function startTheta(net::HybridNetwork,p::Phylip; lbound=0.00001::Float64,factor=2.0::Float64)
    q=extractNQuartets(net,p)
    theta=startTheta(q,net,lbound=lbound,factor=factor) 
    return theta
end

function startTheta(q::Array{Nquartets, 1},net::HybridNetwork; lbound=0.00001::Float64,factor=2.0::Float64)
    
    function findTheta(theta::Float64)
        newT=deepcopy(q)
        newT=pushMomEstlength(newT,theta)#push momest branch lengths for each quartet
        newT=pushMomBL(newT,net)#push averaged momest branch length for the whole tree
        newT=getLikelihoodTh(newT,theta,4/3)#Then use theta, mombl, gamma; get likelihoods for each quartet
        Pseudolikelihood=getLogPseudoLik(newT)#Then compute pseudolikelihood
        return Pseudolikelihood
    end

    lower,upper=upperLower(q,lbound=lbound,factor=factor)
    lower1,upper1=goldenSectionSearch(findTheta,lower,upper,0.01)
    res=Optim.optimize(findTheta,lower1,upper1)
    theta=res.minimizer
    
    return theta
end