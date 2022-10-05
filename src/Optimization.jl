#written by Sungsik Kong 2021-2022

# 6.Optimization

"""
    parentChild

Function to get a vector of the parent-child relationships of node numbers using the node numbers assigned 
in HybridNetwork then transform it into our numbering system (e.g., node number -2 represents the root in 
HybridNetwork but we use 1). The length of the vector must be equal to the number of speciation times in the
topology (or the number of non-leaf vertices) that equals to J=n+h-1 where n represents the number of leaves
and h represents the number of reticulations.

### Input
**`net`**       A tree/network in Type object PhyloNetworks.HybridNetwork\\
"""
function parentChild(net::HybridNetwork) 
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


#function parentChild1(net::HybridNetwork) 
#    parentchild=Array[]
#    for n in net.node
#        while n.number<0
#            parent=parentnode(n.number,net)
#            if parent===nothing
#                push!(parentchild,[0,n.number])
#                @debug "Node $(n.number) is a root."
#            else
#                push!(parentchild,[parent,n.number])
#                @debug "Node $(n.number) is not a root."
#            end
#            break
#        end
#    end
#    #organize array
#    parentchild=sort(parentchild, rev=true)
#    @debug "ParentChild array successfully sorted."
#    #change that negative tau number of HybridNetwork to our format using tauNum
#    for element in parentchild
#        i=element[1]
#        j=element[2]
#        if i!==0
#            element[1]=tauNum(i)
#            @debug "Node number $i changed to $(element[1]) using tauNum."
#        end
#        if j!==0
#            element[2]=tauNum(j)
#            @debug "Node number $j changed to $(element[2]) using tauNum."
#        end
#    end
#    return parentchild
#end

"""
    initialTaus

Using the averaged speciation times and the parent-child relationships, this function ultimately transforms speciation times
for unconstraint optimization of parameters. 

"""
function initialTaus(ave::Array, parentchild::Array)
    
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

#Get starting points for gamma and the edge that has that gamma
"""
    extractGammas

Given the HybridNetwork, this function returns all reticulation edges and its associted gamma.
"""
function extractGammas(net::HybridNetwork)
    obsgam=Float64[]
    retedge=[]
    
    if net.numHybrids==0
        push!(obsgam,1.0)
    else
        for e in net.edge
            if(e.hybrid)
                parentnode=GetParent(e).number
                childnode=GetChild(e).number
                if isempty(obsgam)
                    push!(obsgam,e.gamma)
                    push!(retedge,[parentnode,childnode])
                else
                    redun=false
                    for i in 1:length(retedge)
                        if retedge[i][2]==childnode
                            redun=true
                        end 
                    end
                    if(!redun)
                        push!(obsgam,e.gamma)
                        push!(retedge,[parentnode,childnode])
                        break
                    end
                end
            else
                e.gamma==1.0 || error("Tree edge, but with gamma of $(e.gamma)")
            end
        end
    end

    return obsgam,retedge
end

#get PTGamma *not to be used for optimization
"""
    gamArray

Given the HybridNetwork, it returns the gammas for each parental tree.
"""
function gamArray(net::HybridNetwork)
    obsgam=Float64[]
    retedge=[]
    if net.numHybrids==0
        push!(obsgam,1.0)
    else
        for e in net.edge
            if(e.hybrid)
                parentnode=GetParent(e).number
                childnode=GetChild(e).number
                if isempty(obsgam)
                    push!(obsgam,e.gamma)
                    push!(retedge,[parentnode,childnode])
                else
                    for i in 1:length(retedge)
                        if retedge[i][2]==childnode
                            push!(obsgam,(1-obsgam[i]))
                            push!(retedge,[parentnode,childnode])
                            break
                        else
                            push!(obsgam,e.gamma)
                            push!(retedge,[parentnode,childnode])
                            break
                        end
                    end
                end
            else
                e.gamma==1.0 || error("Tree edge, but with gamma of $(e.gamma)")
            end
        end
    end

    #obsgam,retedge #gamma1,gamma2...reticulation edge with gamma1, edge with gamma2...

    gammas=Float64[]
    dispT=displayedTrees(net,0.0,multgammas=true)
    for t in dispT
        tgammas=Float64[]
        for e in t.edge
            if e.gamma!==1.0
                parentnode=GetParent(e).number
                #childnode=getChild(e).number
                for i in 1:length(retedge)
                    if retedge[i][1]==parentnode
                        push!(tgammas,obsgam[i])
                    end
                end
            end
        end
        push!(gammas,prod(tgammas))
    end
    
    return gammas
end

"""
    gamArray

A better function than the previous one. This one can be used for optimization process.
Given the HybridNetwork, observed gammas in the networks edges, and the edge information,
this function returns inhertiance probabilities for each parental tree.
"""
#get an array for gamma for each PT *used for optimization
function gamArray(net::HybridNetwork,obsgam::Array,retedge::Array)
    #get updated obsgam [gam1,gam2,1-gam1,1-gam2...] and retedge [[edgewithgam1],[edgewith1-gam1],...]
    for e in net.edge
        parentnode=GetParent(e).number
        childnode=GetChild(e).number
        if [parentnode,childnode] in retedge continue
        else
            for i in 1:length(retedge)
                if retedge[i][2]==childnode
                    push!(obsgam,(1-obsgam[i]))
                    push!(retedge,[parentnode,childnode])
                end
            end
        end
    end
    #now obsgam and retedge contains all reticulation edges and corresponding gammas in retedge and obsgam

    gammas=Float64[]
    t=displayedTrees(net,0.0,multgammas=true)#get parenta trees from net
    for pt in t
        ptgam=Float64[]
        for pte in pt.edge
            if pte.gamma!==1.0
                parentnode=GetParent(pte).number
                for i in 1:length(retedge)
                    if retedge[i][1]==parentnode
                        push!(ptgam,obsgam[i])
                    end
                end
            end
        end
        push!(gammas,prod(ptgam))
    end
    
    return gammas
end

function getPseudolik(Optq::Array{Nquartets, 1},theta1::Float64,alpha::Float64,gammas::Array)

    pseudolik=BigFloat[]
    for n in 1:length(Optq[1].nquartet)  #length(q[1].nquartet)=numQuartets
        pnet=Float64[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        for m in 1:length(Optq)
            t1=Optq[m].nquartet[n].mombl[1][1]
            t2=Optq[m].nquartet[n].mombl[1][2]
            t3=Optq[m].nquartet[n].mombl[1][3]

            if Optq[m].nquartet[n].symtype[1]==0
                p=gammas[m]*GetTrueProbsSymm(t1,t2,t3,theta1,alpha)
                pnet=pnet+p
            elseif Optq[m].nquartet[n].symtype[1][1]==1
                p=gammas[m]*GetTrueProbsNetTypes(1,t1,t2,t3,theta1,alpha)
                pnet=pnet+p
            elseif Optq[m].nquartet[n].symtype[1][1]==2
                p=gammas[m]*GetTrueProbsNetTypes(2,t1,t2,t3,theta1,alpha)
                pnet=pnet+p
            elseif Optq[m].nquartet[n].symtype[1][1]==3
                p=gammas[m]*GetTrueProbsNetTypes(3,t1,t2,t3,theta1,alpha)
                pnet=pnet+p
            elseif Optq[m].nquartet[n].symtype[1][1]==4
                p=gammas[m]*GetTrueProbsNetTypes(4,t1,t2,t3,theta1,alpha)
                pnet=pnet+p
            end
        end
                    
        weights=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
        myprobvec1=pnet
        probs = (weights .* myprobvec1)    
        sitefreq=Optq[1].nquartet[n].mspcountsNET[1]
        
        if (sum(ismissing.(probs))!=0) 
            append!(pseudolik,0)
        elseif (sum(any(t->t<=0,probs))==0)
            lik=BigFloat(1.0)
            for i in 1:15
                lik*=BigFloat(probs[i])^sitefreq[i]
            end          
            append!(pseudolik,lik)
        else 
            append!(pseudolik,0)
        end 
    end
    return prod(pseudolik)
end

"""
    Optimization
    
    hehe
"""
Optimization(net::HybridNetwork,p::Phylip)=Optimization(net,p,1000,true)
Optimization(net::HybridNetwork,p::Phylip,NumIter::Int)=Optimization(net,p,NumIter,true)
function Optimization(net::HybridNetwork,p::Phylip,NumIter::Int,printest::Bool;lbound=0.00001,factor=2.0)
    
    q=extractNQuartets(net,p)
    
    #get params
    params=Float64[]
    theta=startTheta(q,net,lbound=lbound,factor=factor) 
    pushMomEstlength(q,theta)
    ave=getmomBLAve(q,net)
    pc=parentChild(net)
    initTau=initialTaus(ave,pc)
    for element in initTau
        push!(params,element)
    end
    
    #get and push gamma to params
    obsgam,retedge=extractGammas(net) #get r1,r2,r3...
    for element in obsgam
        push!(params,asin(sqrt(element)))
    end

    #push theta to params
    push!(params,log(theta))
    
    function objective(params::Array)

        Optq=deepcopy(q) #copy nearly complete Nquartet
        Taus=backTransform(ave,params,pc)
        theta1=exp(params[length(params)])  
        if theta1<1e-5
            theta1=0.02
        end
        if net.numHybrids==0
            gammas=[1.0,0.0]
        else
            gam=Float64[]
            #println(obsgam,retedge)
            for n in 1:net.numHybrids
                push!(gam,sin(params[length(params)-n])^2)
            end
            x,retedges=extractGammas(net) #get r1,r2,r3...
            #println(obsgam,retedge)
            gammas=gamArray(net,gam,retedges)#get gam for pts
        end        
        #println(Taus)
        #if Taus[1]>0
        Optq=pushMomBL(Optq,Taus)
        lq=getPseudolik(Optq,theta1,4/3,gammas)
        y=-1.0*(log(lq))
        
        return y
        #else
        #    return Inf
        #end

    end

    res=Optim.optimize(objective,params,BFGS(linesearch=LineSearches.BackTracking()),Optim.Options(iterations = NumIter))

    params=res.minimizer
    theta1=exp(params[length(params)])  
    theta2=theta1*2
    Taus=backTransform(ave,params,pc)
    if net.numHybrids==0
        gammas=[1.0,0.0]
    else
        obsgam=Float64[]
        for n in 1:net.numHybrids
            push!(obsgam,sin(params[length(params)-n])^2)
        end
        #println("obsgam:$obsgam")
        gammas=obsgam
    end        
        if(printest)        
        println("Optimized Parameter set:\nTaus=$Taus\nGammas=$obsgam\nTheta=$theta2")
    end
 
    return res.minimum,Taus,gammas,theta2,retedge

end