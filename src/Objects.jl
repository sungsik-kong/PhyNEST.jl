#written by Sungsik Kong 2021-2022

#=Set the precision (in bits) to be used for T arithmetic. See https://docs.julialang.org/en/v1/base/numbers/#Base.MPFR.setprecision=#
setprecision(64)

"""
    INPUT

An abstract type that cannot be instantiated, and serve only as nodes in the type graph, with descendant type Phylip.
"""
abstract type INPUT end

"""
    Phylip

Subtype of abstract `INPUT` type with the following attributes:

- `filename`      Name of the input phylip alignment file\\
- `time`          Time taken for parsing the input in seconds\\
- `numtaxa`       Number of taxa given at the first line of the input file\\
- `seqleng`       Sequence length given at the first line of the input file\\
- `nametaxa`      Sequence names given in the input file\\
- `counttaxa`     Unique integer identifier given to each individual in the order of appearance in the input file\\
- `allquartet`    All combinations of quartets for counttaxa\\
- `index`         Just convert allquartet elements into arbitrary index numbers\\
- `spcounts`      Arrays of 15 site pattern frequencies for each quartet in allquartet
"""
mutable struct Phylip <: INPUT
    filename::String 
    time::Float64 
    numtaxa::Integer
    seqleng::Integer
    nametaxa::Array{String,1}
    counttaxa::Array{Integer,1}
    allquartet::Array{Vector{Integer},1}
    index::Array{Vector{String},1}
    spcounts::Array{Vector{Float64},1}
    Phylip()=new("",0.0,0,0,[],[],[],[],[])
    Phylip(inputfile::AbstractString)=new(inputfile,0.0,0,0,[],[],[],[],[])
end

#Now we can see a little nice summary of Phylip in a few lines instead of scary long long waterfall of green lines.
function Base.show(io::IO, p::Phylip)
    disp = "Summary of Phylip File\n"
    disp = disp*"Parsing the file [$(p.filename)] took $(p.time) seconds. \n"
    disp = disp*"Number of taxa: $(p.numtaxa) \n"
    disp = disp*"Species names: $(p.nametaxa) \n"
    disp = disp*"Alignment length (b.p): $(p.seqleng) \n"
    disp = disp*"Site patterns frequencies for $(length(p.spcounts)) quartets computed and stored. \n"
    disp = disp*"Try `show_sp()` function to see all quartet site patterns."
    println(io, disp)
end

"""
    Quartet

An abstract type that cannot be instantiated, and serve only as nodes in the type graph, with descendant type quartets and Network.
"""
abstract type Quartet end

"""
quartets

Subtype of abstract `Quartet` type with the following attributes:

`number`        List of individuals quartets in the order of i,j,k,l\\
`displayed_tree`nth displayed tree that the quartet was extracted from\\
`quartet`       List of quartets in the order of i,j,k,l using the leaf numbers in HybridNetwork\\
`tquartet`      List of quartets in the order of i,j,k,l using the leaf numbers in Phylip\\
`gamma`         Inhertiance probability information provided in HybridNetwork\\
`mspcountsNET`  We can directly use this counts for likelihood calculation.\\
`mrca`          List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\\
`ntau`          Unique number of the taus in the tau used in branchlengths\\
`momestlength`  Branch length identified for each quartet given the mspcounts using moment estimator\\
`average_mom_est_bl`         Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\\
`symtype`       Type of each quartet. It can be either symmetric (type 0) or asymmetric. Asymmetric quartets have four possible topologies:
- Type 1: (i,((j,k),l));
- Type 2: ((i,(j,k)),l); 
- Type 3: (i,(j,(k,l))); 
- Type 4: (((i,j),k),l).\\
`logLik`        It's just a negative likelihood for the quartet in interest
"""
mutable struct quartets <: Quartet
    number::Int64
    displayed_tree::Int64
    quartet::NTuple{}
    tquartet::NTuple{}
    gamma::Float64 
    mspcountsNET::Vector{Float64}
    mrca::NTuple{}
    ntau::NTuple{} 
    momestlength::NTuple{}
    average_mom_est_bl::NTuple{} #used when optimizing.
    symtype::Int64
    logLik::BigFloat #-log...change to just likelihood
    quartets()=new(0,0,(),(),0.0,[],(),(),(),(),9,0.0)
    quartets(number,displayed_tree)=new(number,displayed_tree,(),(),0.0,[],(),(),(),(),9,0.0)
end

"""
Network

Subtype of abstract `Quartet` type with the following attributes:

- `leafname`    Name of the leaves in the order of appearance in the input newick\\
- `leafnumber`  Assigned number for each leaf when reading in the input newick\\
- `gamma`       Gamma as appears on the newick, although irrelevant for our computation\\
- `theta`       Estimated theta using the site pattern frequencies. If Network object is\\
                created without Phylip, theta is not estimated byt default theta=0.001 is displayed.\\
- `quartet`     An array of all quartets in the given topology
"""
mutable struct Network <: Quartet
    leafname::Array{String,1}
    leafnumber::Array{Int64,1}
    gamma::Array{Float64,1}
    theta::Float64
    quartet::Array{quartets,1}
    Network()=new([],[],[],0.001,[])
    Network(leafname,leafnumber)=new(leafname,leafnumber,[],0.001,[])
end

#Again, we don't want to show all the scary attributes but only a summary of the extracted quartets.
function Base.show(io::IO, nq::Network)
    disp = "Summary of parsed $(typeof(nq))\n"
    disp = disp*"$(length(nq.gamma)) displayed tree(s) present in the network.\n"
    disp = disp*"Inheritance parameter for each displayed tree: $(nq.gamma) \n"
    #disp = disp*"Estimated theta: $(nq.theta) \n"
    disp = disp*"Leaf name: $(nq.leafname) \n"
    disp = disp*"Leaf number: $(nq.leafnumber) \n"
    disp = disp*"Total number of quartets: $(length(nq.quartet))"
    println(io, disp)
end
