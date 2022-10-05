#written by Sungsik Kong 2021-2022

# 0.Preliminaries

setprecision(64)


#Set the precision (in bits) to be used for T arithmetic.
#see https://docs.julialang.org/en/v1/base/numbers/#Base.MPFR.setprecision

abstract type INPUT end
"""
    Phylip

Subtype of abstract `INPUT` type. Sequence alignment with the following attributes:

- `filename`    The actual filename of the input sequence alignment\\
- `time`    Time it took to intake the filename in seconds\\
- `numtaxa` The number of taxa appears on the first line of phylip\\
- `seqleng` The sequence length appears on the first line of phylip\\
- `nametaxa`    Names of each individual in the order of appearance in phylip\\
- `counttaxa`   Numbering each individual in the order of appearance in phylip\\
- `allquartet`  All combinations of quartets using nametaxa\\
- `index`   Just convert allquartet elements into arbitrary index numbers\\
- `spcounts`    Arrays of 15 site pattern frequencies 
"""
mutable struct Phylip <: INPUT
    filename::String #the actual filename of the input sequence alignment
    time::Float64 #time it took to intake the $filename in seconds
    numtaxa::Integer #the number of taxa appears on the first line of phylip
    seqleng::Integer #the sequence length appears on the first line of phylip
    nametaxa::Array{String,1} #names of each individual in the order of appearance in phylip
    counttaxa::Array{Integer,1} #numbering each individual in the order of appearance in phylip
    allquartet::Array{Vector{Integer},1} #all combinations of quartets using nametaxa
    index::Array{Vector{String},1} #just convert allquartet elements into arbitrary index numbers
    spcounts::Array{Vector{Float64},1} #arrays of 15 site pattern frequencies 
    #the message pushed to the filename attribute in the object Phylip when no input filename was specified.
    #Phylip()=@error "No input sequence file specified."
    Phylip()=new("",0.0,0,0,[],[],[],[],[])
    Phylip(input::AbstractString)=new(input,0.0,0,0,[],[],[],[],[])
end

#This will allow *not* showing the scary long elements of the object **Phylip**, instead it will summarize some relevant information.
function Base.show(io::IO, p::Phylip)
    disp = "Summary: $(typeof(p)) File \n"
    disp = disp*"Parsing the file [$(p.filename)] took [$(p.time)] seconds. \n"
    disp = disp*"Number of taxa: $(p.numtaxa) \n"
    disp = disp*"Species names: $(p.nametaxa) \n"
    disp = disp*"Sequence length (b.p): $(p.seqleng) \n"
    disp = disp*"Site patterns frequencies for [$(length(p.spcounts))] quartets computed and stored."
    println(io, disp)
end


abstract type NQuartet end
"""
    nquartets

Subtype of abstract `NQuartet` type. Stores relevant information of a quartet extracted from a topology with the following attributes:

- `number`  List of individuals quartets in the order of i,j,k,l\\
- `tquartet`    List of quartets in the order of i,j,k,l\\
- `indexNET`    Index number for each quartet in the order appears in quartnet. \\
- `mspcountsNET`    We can directly use this counts for likelihood calculation.\\
- `mrca`    List of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl\\
- `ntau`    Unique number of the taus in the tau used in branchlengths\\
- `branchlength`    Branch length from the mrca to tau=0 in the order of i and j (ij),ik,il,jk,jl,and kl\\
- `momestlength`    Branch length identified for each quartet given the mspcounts using moment estimator - Will get filled later because theta is required\\
- `mombl`   Branch length that is averaged for the entire tree/network - Will get filled later because theta is required\\
- `symtype` Type of each quartet. It can be either symmetric or asymmetric. Asymmetric (type 0) can have four possible topologies and we define type 1 as (i,((j,k),l)); type 2 as ((i,(j,k)),l); type 3 as (i,(j,(k,l))); and type 4 as (((i,j),k),l).\\
- `logLik`  It's just a likelihood for the quartet
"""
mutable struct nquartets <: NQuartet
    number::Int64
    quartet::Array{Array{Int64,1}} #list of individuals quartets in the order of i,j,k,l
    tquartet::Array{Array{Int64,1}} #list of quartets in the order of i,j,k,l
    indexNET::Array{Vector{String},1} #index number for each quartet in the order appears in quartnet. 
    mspcountsNET::Array{Vector{Float64},1} #We can directly use this counts for likelihood calculation.
    mrca::Array{Array{Int64,1}} #list of common ancesters of two taxa in the order of i and j (ij),ik,il,jk,jl,and kl
    ntau::Array{} #unique number of the taus in the tau used in branchlengths
    branchlength::Array{Array{Float64,1}} #branch length from the mrca to tau=0 in the order of i and j (ij),ik,il,jk,jl,and kl
    momestlength::Array{} #branch length identified for each quartet given the mspcounts using moment estimator *Will get filled later because theta is required
    mombl::Array{} #branch length that is averaged for the entire tree/network *Will get filled later because theta is required
    symtype::Array{Int64,1} #type of each quartet. It can be either symmetric or asymmetric. Asymmetric (type 0) can have four
        #possible topologies and we define type 1 as (i,((j,k),l)); type 2 as ((i,(j,k)),l); type 3 as (i,(j,(k,l))); and 
        #type 4 as (((i,j),k),l).
    #pquartet::Array{Array{Int64,1}} #Transform the members in .quartet attribute to match with the phylip numbering %We don't need it anymore
    #index::Array{Vector{String},1} #index number for each quartet so we can find the site pattern counts accordingly. This index number will
        #not be i*j*k*l for the asymmetric tree types 1,2,and 4. For these types, we will have to use i*l*j*k for type 1,
        #l*i*j*k for type 2, and l*k*i*j for type 4. %We don't need it anymore
    #mspcounts::Array{Vector{Float64},1} %We don't need it anymore
    logLik::BigFloat #-log
    
    nquartets(number)=new(number,[],[],[],[],[],[],[],[],[],[],0.0)
    #nquartets(number)=new(number,[],[],[],[],[],[],[],[],[],[],[],[],0.0)
end


"""
    Nquartets

Subtype of abstract `NQuartet` type. Stores relevant information of a tree extracted from a topology with the following attributes:

- `leafname`    Name of the leaves as appears on the newick\\
- `leafnumber`  Assigned number for each leaf when reading in the newick\\
- `gamma`   Gamma as appears on the newick, although irrelevant for our computation\\
- `nquartet`    An array of all quartets in the given topology
"""
mutable struct Nquartets <: NQuartet
    leafname::Array{String,1}
    leafnumber::Array{Int64,1}
    gamma::Float64
    nquartet::Array{nquartets,1}
    Nquartets(leafname,leafnumber)=new(leafname,leafnumber,0.5,[])
end

#Again, we don't want to show all the scary attributes so we only display a summary of extracted quartets and some relevant information for a topology. 
function Base.show(io::IO, nq::Nquartets)
    disp = "Summary of parsed $(typeof(nq))\n"
    disp = disp*"Gamma: $(nq.gamma) \n"
    disp = disp*"Leafname: $(nq.leafname) \n"
    disp = disp*"Number of Quartets per tree: $(length(nq.nquartet))"
    println(io, disp)
end

