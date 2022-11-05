#testing functions in src/Quartets.jl

@testset "Testing quartet extraction" begin
    p=readPhylipFile!("../example/quartet_500.txt")
    @test PhyNEST.tauNum(-2)==1
    #A quartet tree case
    #symmetric quartet
    x=readTopology("((1,2),(3,4));") 
    leafname,leafnumber=PhyNEST.getLeafInfo(x)
    @test leafname==["1", "2", "3", "4"]
    @test leafnumber==[1,2,3,4]
    q=Nquartets(leafname,leafnumber)#creates Nquartets type.
    nq=PhyNEST.listAllQs(q)
    @test nq.nquartet[1].quartet==[[1,2,3,4]]
    @test PhyNEST.childnode(-2,x)==-3
    @test PhyNEST.childnode(-3,x)==1
    @test PhyNEST.childnode(-4,x)==3
    for i in x.node if(i.leaf) @test PhyNEST.childnode(i.number,x)===nothing end end
    @test PhyNEST.parentnode(1,x)==-3
    @test PhyNEST.parentnode(2,x)==-3
    @test PhyNEST.parentnode(3,x)==-4
    @test PhyNEST.parentnode(4,x)==-4
    @test PhyNEST.parentnode(-4,x)==-2
    @test PhyNEST.parentnode(-3,x)==-2
    @test PhyNEST.parentnode(-2,x)===nothing
    @test PhyNEST.mrca(1,2,x)==-3
    @test PhyNEST.mrca(1,3,x)==-2
    @test PhyNEST.mrca(1,4,x)==-2
    @test PhyNEST.mrca(2,3,x)==-2
    @test PhyNEST.mrca(2,4,x)==-2
    @test PhyNEST.mrca(3,4,x)==-4
    PhyNEST.mrcapush(x,nq.nquartet[1])
    @test nq.nquartet[1].mrca==[[-3, -2, -2, -2, -2, -4]]
    nq=PhyNEST.symType(nq.nquartet[1])
    @test nq.symtype==[0]
    PhyNEST.uniqtaus(nq)
    @test nq.ntau==[[3,2,1]]
    dp=PhyNEST.DictionaryP(p)
    for i in p.nametaxa @test haskey(dp,i) end
    dq=PhyNEST.DictionaryQ(q)
    for i in leafname @test haskey(dq,i) end
    dn=PhyNEST.DictionaryN(q,p)
    for i in leafnumber @test haskey(dn,i) end
    PhyNEST.pquartet(q,dn)
    @test q.nquartet[1].tquartet==[[3, 4, 2, 1]]
    PhyNEST.binaryIndexforQuartet(q)
    @test q.nquartet[1].indexNET==[["00000011000001000000001000000001"]]
    PhyNEST.moveSPcounts(q,p)
    @test length(q.nquartet[1].mspcountsNET[1])==15
    #type 1
    x=readTopology("(1,((2,3),4));")
    leafname,leafnumber=PhyNEST.getLeafInfo(x)
    @test leafname==["1", "2", "3", "4"]
    @test leafnumber==[1,2,3,4]
    q=Nquartets(leafname,leafnumber)#creates Nquartets type.
    nq=PhyNEST.listAllQs(q)
    @test nq.nquartet[1].quartet==[[1,2,3,4]]
    @test PhyNEST.childnode(-2,x)==1
    @test PhyNEST.childnode(-3,x)==-4
    @test PhyNEST.childnode(-4,x)==2
    for i in x.node if(i.leaf) @test PhyNEST.childnode(i.number,x)===nothing end end
    @test PhyNEST.parentnode(1,x)==-2
    @test PhyNEST.parentnode(2,x)==-4
    @test PhyNEST.parentnode(3,x)==-4
    @test PhyNEST.parentnode(4,x)==-3
    @test PhyNEST.parentnode(-4,x)==-3
    @test PhyNEST.parentnode(-3,x)==-2
    @test PhyNEST.parentnode(-2,x)===nothing
    @test PhyNEST.mrca(1,2,x)==-2
    @test PhyNEST.mrca(1,3,x)==-2
    @test PhyNEST.mrca(1,4,x)==-2
    @test PhyNEST.mrca(2,3,x)==-4
    @test PhyNEST.mrca(2,4,x)==-3
    @test PhyNEST.mrca(3,4,x)==-3  
    PhyNEST.mrcapush(x,nq.nquartet[1])
    @test nq.nquartet[1].mrca==[[-2, -2, -2, -4, -3, -3]]
    nq=PhyNEST.symType(nq.nquartet[1])
    @test nq.symtype==[1]
    PhyNEST.uniqtaus(nq)
    @test nq.ntau==[[3,2,1]]
    dp=PhyNEST.DictionaryP(p)
    for i in p.nametaxa @test haskey(dp,i) end
    dq=PhyNEST.DictionaryQ(q)
    for i in leafname @test haskey(dq,i) end
    dn=PhyNEST.DictionaryN(q,p)
    for i in leafnumber @test haskey(dn,i) end
    PhyNEST.pquartet(q,dn)
    @test q.nquartet[1].tquartet==[[3, 4, 2, 1]]
    PhyNEST.binaryIndexforQuartet(q)
    @test q.nquartet[1].indexNET==[["00000011000001000000001000000001"]]
    PhyNEST.moveSPcounts(q,p)
    @test length(q.nquartet[1].mspcountsNET[1])==15
    #single reticulation n5h1 case

    #double reticulation n7h1 case

end

#getLeafInfo
#listAllQs
#childnode
#parentnode
#mrca
#mrcapush
#symType
#uniqtaus
#DictionaryP
#DictionaryQ
#DictionaryN
#pquartet
#binaryIndexforQuartet
#moveSPcounts
####findGamma
#extractNQuartets(net::HybridNetwork)
#extractNQuartets(net::HybridNetwork,p::Phylip)
#printQuarts