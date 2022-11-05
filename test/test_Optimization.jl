#testing functions in src/Optimization.jl

@testset "Testing Optimization" begin
    #Parameterization and Transformations
    x=readTopology("(5,(4,((1,(2)#H6:::0.1),(3,#H6:::0.9))));")
    p=readPhylipFile!("../example/n5h1_5k.txt")
    n=x.numTaxa
    h=x.numHybrids

    ##SPECIATION TIMES #Function parentChild
    parentchild=PhyNEST.parentChild(x)#number of speciation times==n+h-1
        @test length(parentchild)==(n+h-1)
    childs=Integer[]#every child must be unique
    for i in parentchild push!(childs,i[2]) end
        @test length(parentchild)==length(childs)
        @test parentchild[1][1]==0#the first one in the array should be the root
        root=0#there must only one root node
    for i in parentchild if i[1]==0; root+=1 end end
        @test root==1 

    #Function initialTaus
    q=extractNQuartets(x,p)
    theta=PhyNEST.startTheta(q,x,lbound=1e-5,factor=2.0) 
    PhyNEST.pushMomEstlength(q,theta)
    #getmomBLAve see if all averaged speciations times are positives
    ave=PhyNEST.getmomBLAve(q,x)
    for aveTau in ave if !isnan(aveTau) @test aveTau >= 0 end end
    #initialTaus
    initTau=PhyNEST.initialTaus(ave,parentchild)
    @test length(initTau)==(n+h-1)
    #backTransform
    Taus=PhyNEST.backTransform(ave,initTau,parentchild)
    #check if back-transformed speciation times are same as averaged speciation times
    @test length(Taus)==length(ave)
    for i in 1:length(Taus)
        if !isnan(ave[i]) @test isapprox(Taus[i],ave[i],atol=1e-5) end
    end
    
    ##INHERITANCE PROBABILITY
    obsgam,retedge=PhyNEST.extractGammas(x)
    @test length(obsgam)==length(retedge)==h
    for gam in obsgam @test gam>=0 && gam<=1.0 end #gamma must be (0,1)
    @test sum(PhyNEST.gamArray(x))==1.0 #see if the inheritance probabilities for each parental trees adds up to 1
    PTgam=PhyNEST.gamArray(x,obsgam,retedge)
    @test sum(PTgam)==1.0
    @test length(PTgam)==2^h

    #LIKELIHOODS and OPTIMIAZATION
    Optq=PhyNEST.pushMomBL(q,Taus)
    lq=PhyNEST.getPseudolik(Optq,theta,4/3,PTgam)
    @test lq>0.0 && lq<Inf && !isnan(lq) #pseudolik is some number
    Optmin,OptTaus,Optgammas,Opttheta,Optretedge=Optimization(x,p,1000,false)
    @test Optmin>0.0 && Optmin<Inf && !isnan(Optmin) #pseudolik is some number
    for i in OptTaus if i!==0.0 @test i>0.0 end end
    for i in Optgammas @test i>=0.0 && i<=1.0 end
    @test Opttheta>0.0 && Opttheta<1.0
    for i in Optretedge @test i in retedge end
end