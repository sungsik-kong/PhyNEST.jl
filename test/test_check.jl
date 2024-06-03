rep=100#number of reps


@testset "Testing misc.function get_average" begin
    for i in 1:rep
        i=[0,10,50,40,0.5]
        ave=PhyNEST.get_average(i)
        @test ave==20.1
    end
end

#testing functions in src/Probabilities.jl
@testset "Testing site pattern proabilities" begin
    for i in 1:rep
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
        @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsAsymmTypes(1,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
        @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsAsymmTypes(2,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
        @test isapprox(sum(wgt .* PhyNEST.GetTrueProbsAsymmTypes(4,t1,t2,t3,theta,alpha)),1.0,atol=1e-5)
        #sum of simulated site pattern frequencies must equal to specified length
        @test sum(sim_sp_freq(0,t1,t2,t3,theta,alpha,n))==n
        @test sum(sim_sp_freq(1,t1,t2,t3,theta,alpha,n))==n
        @test sum(sim_sp_freq(2,t1,t2,t3,theta,alpha,n))==n
        @test sum(sim_sp_freq(3,t1,t2,t3,theta,alpha,n))==n
        @test sum(sim_sp_freq(4,t1,t2,t3,theta,alpha,n))==n
    end
end

#testing functions in src/MomentEstimat.jl
@testset "Testing Method-of-moment estimates of speciation times" begin
    for i in 1:rep
        #drawing random parameters...
        #constraints: t1<t2<t3; 1e-5<theta<0.1; alpha=4/3
        t1=rand(1:.0001:3)
        t2=rand(4:.0001:6)
        t3=rand(6:.0001:10)
        theta=rand(0.00001:.0001:0.1)
        alpha=4/3
        n=Int(rand(1e5:1e7))
        tol=0.2 #setting this tolerance too low would result in test failures
        spcounts0=sim_sp_freq(0,t1,t2,t3,theta,alpha,n)
        est1,est2,est3=PhyNEST.mom_est_sym(spcounts0,theta)
        @test isapprox(est1,t1,atol=tol)   
        @test isapprox(est2,t2,atol=tol)   
        @test isapprox(est3,t3,atol=tol)   
        spcounts1=sim_sp_freq(1,t1,t2,t3,theta,alpha,n)
        est1,est2,est3=PhyNEST.mom_est_asym(1,spcounts1,theta)
        @test isapprox(est1,t1,atol=tol)   
        @test isapprox(est2,t2,atol=tol)   
        @test isapprox(est3,t3,atol=tol)   
        spcounts2=sim_sp_freq(2,t1,t2,t3,theta,alpha,n)
        est1,est2,est3=PhyNEST.mom_est_asym(2,spcounts2,theta)
        @test isapprox(est1,t1,atol=tol)   
        @test isapprox(est2,t2,atol=tol)   
        @test isapprox(est3,t3,atol=tol)   
        spcounts3=sim_sp_freq(3,t1,t2,t3,theta,alpha,n)
        est1,est2,est3=PhyNEST.mom_est_asym(3,spcounts3,theta)
        @test isapprox(est1,t1,atol=tol)   
        @test isapprox(est2,t2,atol=tol)   
        @test isapprox(est3,t3,atol=tol)   
        spcounts4=sim_sp_freq(4,t1,t2,t3,theta,alpha,n)
        est1,est2,est3=PhyNEST.mom_est_asym(4,spcounts4,theta)
        @test isapprox(est1,t1,atol=tol)   
        @test isapprox(est2,t2,atol=tol)   
        @test isapprox(est3,t3,atol=tol)   
    end
end