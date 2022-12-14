Dstat(outgroup::String,p::Phylip)=Dstat(outgroup::String,p::Phylip,0.05,false)
Dstat(outgroup::String,p::Phylip,pval::Float64)=Dstat(outgroup::String,p::Phylip,pval::Float64,false)
function Dstat(outgroup::String,p::Phylip,pval::Float64,displayall::Bool)
    
    outgroup1=0
    ndist=Normal(0,1)
    quartet=[]
    sitepattern=[]

    for n in 1:length(p.nametaxa)
        
        if p.nametaxa[n]==outgroup
            outgroup1=p.counttaxa[n]
            break
        end
    end

    outgroup1!==0 || error("Cannot find the specified outgroup '$outgroup' in Phylip input file.")

    for n in 1:length(p.allquartet)
        if p.allquartet[n][1]==outgroup1
            push!(quartet,p.allquartet[n])
            push!(sitepattern,[p.spcounts[n][7],p.spcounts[n][9]])
        else
            continue
        end
    end

        for n in 1:length(sitepattern)
            ABAB=sitepattern[n][1]
            ABBA=sitepattern[n][2]
            d=(ABBA-ABAB)/(ABBA+ABAB)
            z=d/(2*sqrt((0.25/(ABBA+ABAB))))
            
            pv=1-cdf(ndist,z)
            if displayall==true
                println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv")
            else
                if pv<=(pval)
                    println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv <= SIGNIFICANT!")
                else
                    println("$([p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]), d=$d, z=$z, p=$pv")
                end
            end

        end

    end


###Need to deal with avg_obs=num_obs/(1*1*1*1)#(self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
function HyDe(outgroup::String,p::Phylip; p_value=0.05::Float64, filter=true::Bool)
    outgroup1=0
    quartet=[]
    sitepattern=[]
    num_obs = p.seqleng
     
    

    for n in 1:length(p.nametaxa)
        if p.nametaxa[n]==outgroup
            outgroup1=p.counttaxa[n]
            break
        end
    end

    outgroup1!==0 || error("Cannot find the specified outgroup '$outgroup' in Phylip input file.")

    for n in 1:length(p.allquartet)
        if p.allquartet[n][1]==outgroup1
            push!(quartet,p.allquartet[n])
            push!(sitepattern,[p.spcounts[n][4],p.spcounts[n][7],p.spcounts[n][9],p.spcounts[n][3],p.spcounts[n][2],p.spcounts[n][6]])
        else
            continue
        end
    end


    for n in 1:length(sitepattern)

        p9 = (sitepattern[n][3]+0.05)/num_obs#ABBA
        p7 = (sitepattern[n][2]+0.05)/num_obs#ABAB
        p4 = (sitepattern[n][1]+0.05)/num_obs#AABB

        avg_obs=num_obs/(1*1*1*1)#(self.outIndex.shape[0] * p1.shape[0] * hyb.shape[0] * p2.shape[0])
        avobs=avg_obs
        
        #z-score
        obs_invp1 = avobs * (p9 - p7)
        obs_invp2 = avobs * (p4 - p7)
        obs_var_invp1 = avobs * p9 * (1 - p9) + avobs * p7 * (1 - p7) + 2 * avobs * p9 * p7
        obs_var_invp2 = avobs * p4 * (1 - p4) + avobs * p7 * (1 - p7) + 2 * avobs * p4 * p7
        obs_cov_invp1_invp2 = -1 * avobs * p9 * p4 + avobs * p9 * p7 + avobs * p7 * p4 + avobs * p7 * (1 - p7)
        ratio = obs_invp2 / obs_invp1;
        GH_ts = ((obs_invp1) * (ratio) / sqrt(obs_var_invp1 * (ratio^2) - 2.0 * obs_cov_invp1_invp2 * ratio + obs_var_invp2))

        temp = -99999.9
        if p7 > p9 && p7 < p4
            GH_ts=temp
        elseif GH_ts > -99999.9 && GH_ts < 99999.9
            GH_ts=GH_ts
        else
            GH_ts=temp
        end
        my_z=GH_ts


        #gamma
        _c_num   = avg_obs * (p9 - p7)
        _c_denom = avg_obs * (p4 - p7)
        _c       = _c_num / _c_denom
        gamma = _c / (1 + _c)


        #p-val
        a1 =  0.254829592;
        a2 = -0.284496736;
        a3 =  1.421413741;
        a4 = -1.453152027;
        a5 =  1.061405429;
        px  =  0.3275911;
        sign = 1;
        if my_z < 0 sign = -1 end
        z = abs(my_z) / sqrt(2.0);
        t = 1.0 / (1.0 + px * z);
        y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-z * z);

        p_val=1.0 - (0.5 * (1.0 + sign * y))


        newQ=[p.nametaxa[quartet[n][1]],p.nametaxa[quartet[n][2]],p.nametaxa[quartet[n][3]],p.nametaxa[quartet[n][4]]]
        if filter
            if p_val<(p_value/3)
                println("$newQ, gamma: $gamma, z-score: $GH_ts, p-val: $p_val")
            end
        else
            println("$newQ, gamma: $gamma, z-score: $GH_ts, p-val: $p_val")
        end
    end

end



function LRT(t1,t2,p)
    optt1,tau,gamma=Optimization(t1,p,1000)
    optt2,tau,gamma=Optimization(t2,p,1000)

    supt1=optt1.minimum
    supt2=optt2.minimum

    lr= -2 * log(supt1/supt2)

    return lr 

end



