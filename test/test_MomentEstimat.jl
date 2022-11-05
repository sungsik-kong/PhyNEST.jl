#testing functions in src/MomentEstimat.jl

@testset "Testing Method-of-moment estimates of speciation times" begin
    j=10#number of repeats to test function momentEstimat starting at Ln57
    #drawing random parameters...
    #constraints: t1<t2<t3; 1e-5<theta<0.1; alpha=4/3
    t1=rand(1:.0001:3)
    t2=rand(4:.0001:6)
    t3=rand(6:.0001:10)
    theta=rand(0.00001:.0001:0.1)
    alpha=4/3
    n=Int(rand(1e5:1e7))
    spcounts0=simspcounts(0,t1,t2,t3,theta,alpha,n)
    p1,p2=PhyNEST.spProbSym(spcounts0)
    @test sum(p1)>0
    @test sum(p2)>0
    est1,est2,est3=PhyNEST.momEstSym(spcounts0,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    spcounts1=simspcounts(1,t1,t2,t3,theta,alpha,n)
    p=PhyNEST.MOMspProbAsym1(spcounts1)
    @test sum(p)>0
    est1,est2,est3=PhyNEST.MOMEstAsym(1,spcounts1,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    spcounts2=simspcounts(2,t1,t2,t3,theta,alpha,n)
    p=PhyNEST.MOMspProbAsym2(spcounts2)
    @test sum(p)>0
    est1,est2,est3=PhyNEST.MOMEstAsym(2,spcounts2,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    spcounts3=simspcounts(3,t1,t2,t3,theta,alpha,n)
    p=PhyNEST.MOMspProbAsym3(spcounts3)
    @test sum(p)>0
    est1,est2,est3=PhyNEST.MOMEstAsym(3,spcounts3,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    spcounts4=simspcounts(4,t1,t2,t3,theta,alpha,n)
    p=PhyNEST.MOMspProbAsym4(spcounts4)
    @test sum(p)>0
    est1,est2,est3=PhyNEST.MOMEstAsym(4,spcounts4,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    type=rand(0:4)
    spcounts=simspcounts(type,t1,t2,t3,theta,alpha,n)
    est1,est2,est3=PhyNEST.momentEstimat(type,spcounts,theta)
    @test isapprox(est1,t1,atol=0.1)   
    @test isapprox(est2,t2,atol=0.1)   
    @test isapprox(est3,t3,atol=0.1)   
    

end