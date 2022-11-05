#testing functions in src/Probabilities.jl

@testset "Testing site pattern proabilities" begin
    #drawing random parameters...
    #constraints: t1<t2<t3; 1e-5<theta<0.1; alpha=4/3
    t1=rand(1:.0001:3)
    t2=rand(4:.0001:6)
    t3=rand(6:.0001:10)
    theta=rand(0.00001:.0001:0.1)
    alpha=4/3
    n=Int(rand(1e5:1e7))
    wgt=[4,12,12,12,24,12,12,24,12,12,24,24,24,24,24]
    #sum of probabilities must be approx(1), tolerance set to 1e-5
    @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsSymm(t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
    @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsAsymm(t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
    @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsNetTypes(1,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
    @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsNetTypes(2,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
    @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsNetTypes(4,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
    #sum of simulated site pattern frequencies must equal to specified length
    @test sum(simspcounts(0,t1,t2,t3,theta,alpha,n))==n
    @test sum(simspcounts(1,t1,t2,t3,theta,alpha,n))==n
    @test sum(simspcounts(2,t1,t2,t3,theta,alpha,n))==n
    @test sum(simspcounts(3,t1,t2,t3,theta,alpha,n))==n
    @test sum(simspcounts(4,t1,t2,t3,theta,alpha,n))==n
end
