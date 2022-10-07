#written by Sungsik Kong 2021-2022

# 1.Input

"""
    readPhylipFile!(inputfile::AbstractString, writecsv::Bool)
    readPhylipFile!(inputfile::AbstractString)

Takes in, read, and parse the input phylip file. File name is specified as a string. If `writecsv=true`, list of 
site pattern frequencies from the input phylip for all permutations of four taxa is saved as a `.csv` file in the 
working directory. The `writecsv` boolean is false by default. The function will perform assuming `writecsv=false` 
if not specified.

### Input
`inputfile` Name of phylip file as a AbstractString\\
`writecsv` Boolean to write site pattern frequencies in CSV file\\

"""
function readPhylipFile!(inputfile::AbstractString;writecsv=false::Bool,csvname=""::AbstractString,showProgress=true::Bool)
    try
        p=@timed readPhylipFile(inputfile,writecsv,csvname,showProgress)
        p[1].time=round(p[2],digits=3)
        printCheckPoint(p[1])
        return p[1]
    catch(err)
        display(err)
    end
end
#readPhylipFile!(inputfile::AbstractString;showProgress::Bool)=readPhylipFile!(inputfile;writecsv=false,showProgress=showProgress)
#readPhylipFile!(inputfile::AbstractString;writecsv::Bool)=readPhylipFile!(inputfile;writecsv=writecsv,showProgress=false)

readPhylipFile(inputfile::AbstractString)=readPhylipFile(inputfile,false,"",true)
function readPhylipFile(inputfile::AbstractString,writecsv::Bool,csvname::AbstractString,showProgress::Bool)

    p=Phylip(inputfile)
    UniqueBase,BaseCounts=PhylipFileInfo(inputfile, p, showProgress)
    getUniqueQuartets(p)
    sitePatternCounts(p,UniqueBase,BaseCounts)
    spRearrange(p)
    binaryIndexforQuartet(p)
    
    #write site pattern counts for all quartets into .csv
    if(writecsv) writeSitePatternCounts(p,writecsv,csvname,inputfile) end
    
    return p
end

PhylipFileInfo(inputfile::AbstractString, p::Phylip)=PhylipFileInfo(inputfile,p,true)
function PhylipFileInfo(inputfile::AbstractString, p::Phylip, showProgress::Bool)
    #taxaMatch=false
    seq=String[]
    countTaxa1=0#the number of taxa indicated in the first line
    countTaxa2=0#the number of taxa identified by counting the number of lines
    seqLength=0
    
    #println("\nReport for input Phylip file [$inputfile]:")
    open(inputfile) do following
        line = 0
        for l in eachline(following)
            line += 1 #keep count the number of lines
            l=lstrip(l)#remove white spaces before any Char
            eachword=split(l, isspace)#if there is a space between two columns
            for c in eachword
                if c==""
                    deleteat!(eachword, eachword .== c);
                end
            end
            col=length(eachword)
            #if col!==2#if there is no space between two columns but delimited by a tab
            #     eachword=split(l, "\t")#but the space is caused by neither tab or space, then there will be an error <- honestly this is preferred. Should be in this format.
            #end
            col==2 || error("Cannot parse phylip file [$inputfile] on line $line. Every delimitation must be done by a single space or a tab.")

            if line==1 #for the first line of phylip
                #println("\tThere are $(eachword[1]) taxa with Sequence length of $(eachword[2]).")
                try countTaxa1=parse(Int64, eachword[1])#Transform the string value into an integer
                catch e error("Check the first line of [$inputfile]. Number of individual must be an integer, but received '$(eachword[1])'")
                    break
                end
                
                try seqLength=parse(Int64, eachword[2])#Transform the string value into an integer
                catch e error("Check the first line of [$inputfile]. Sequence length must be an integer, but received '$(eachword[2])'") 
                    break
                end
            else
                #println("\tDNA Sequence of the taxon number $(line-1) with name [$(eachword[1])] is composed of nucleotides $(unique(eachword[2])).")
                countTaxa2+=1
                push!(p.nametaxa,eachword[1])#pushed sequence name to the Type Phylip.
                push!(p.counttaxa,countTaxa2)#we designate a unique integer starting from 1 for each individual
                seqLength==length(eachword[2]) || error("The sequence length indicated in the first line of the file [$inputfile] does not match with actual length for the taxon number $countTaxa2 with name [$(eachword[1])].")
                push!(seq,eachword[2])
            end
        end

        countTaxa1!==0 || error("There is no sequence taht we can parse.")
        countTaxa1==countTaxa2 || error("The number of taxa indicated in the first line of the file [$inputfile] is [$countTaxa1], but there are [$(line-1)] sequences.")
            p.numtaxa=countTaxa1 #Begin to fill in the numtaxa element in Type Phylip
            p.seqleng=seqLength
        
        ###ppbaseUnique
        ppbase=Array[]
        numInd=p.numtaxa

        for eachseq in 1:seqLength
            oneSite=Char[]
            for eachind in 1:numInd
                append!(oneSite,seq[eachind][eachseq])
            end
            push!(ppbase,oneSite)
        end
        UniqueBase=unique(ppbase)
        
        #=using package ProgressMeter
        #BaseCounts=[count(==(element),ppbase) for element in UniqueBase] #original line without using ProgressMeter
        #p = Progress(n, 1, "Parsing [$inputfile]...", 50)#using package ProgressMeter=#
        BaseCounts=Int[]
        
        n = length(UniqueBase)
        if(showProgress)
            p = Progress(n, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
            for f in 1:length(UniqueBase)
                BaseCounts=append!(BaseCounts,count(==(UniqueBase[f]),ppbase))
                next!(p)#using package ProgressMeter
            end
        else
            for f in 1:length(UniqueBase)
                BaseCounts=append!(BaseCounts,count(==(UniqueBase[f]),ppbase))
            end
        end
        
        length(UniqueBase)==length(BaseCounts) || error("Something went wrong while reading the sequences. The length of UniqueBase and BaseCounts are not the same.")
        return UniqueBase, BaseCounts
    end
end

function getUniqueQuartets(p::Phylip)
    
    numind=p.numtaxa
    numind!==0 || error("There is no sequence taht we can parse.")
    typeof(numind)==Int64 || error("Error while getting all posible quartets using Type Phylip . Number of individuals stored is not Integer but $(typeof(numind))")
    taxa=p.counttaxa #list of integers assigned for each taxa stored in Type Phylip
    
    for i in 1:numind
        for j in i+1:numind
            for k in j+1:numind
                for l in k+1:numind
                    push!(p.allquartet, [taxa[i],taxa[j],taxa[k],taxa[l]]); 
                end
            end
        end
    end

    length(p.allquartet)==binomial(numind,4) || error("There are only $(length(p.allquartet)) quartets when there should be $(binomial(numind,4)).")    
end

function sitePatternCounts(p::Phylip,ppbase::Array,counts::Array)
    Allquartet=p.allquartet
    numquartet=length(p.allquartet) 
    numquartet!==0 || error("No quartet found to make site pattern counts.")

    for eachquartet in 1:numquartet
        spcount=Float64[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]  
        i=Allquartet[eachquartet][1]
        j=Allquartet[eachquartet][2]
        k=Allquartet[eachquartet][3]
        l=Allquartet[eachquartet][4]

        #assuming length(counts)==length(p.allquartet)==length(ppbase)...
        ##which is checked in PhylipFileInfo while parsing the sequences.
        for n in 1:length(counts)
            if ppbase[n][i]==ppbase[n][j]==ppbase[n][k]==ppbase[n][l] #AAAA
                spcount[1]=spcount[1]+counts[n]
            elseif ppbase[n][i]==ppbase[n][j]==ppbase[n][k]!==ppbase[n][l] #AAAB
                spcount[2]=spcount[2]+counts[n]
            elseif ppbase[n][i]==ppbase[n][j]==ppbase[n][l]!==ppbase[n][k] #AABA
                spcount[3]=spcount[3]+counts[n]
            elseif ppbase[n][i]==ppbase[n][j]!==ppbase[n][k]==ppbase[n][l] #AABB
                spcount[4]=spcount[4]+counts[n]
            elseif ppbase[n][i]==ppbase[n][j]!==ppbase[n][k]!==ppbase[n][l] #AABC
                spcount[5]=spcount[5]+counts[n]
            elseif ppbase[n][i]==ppbase[n][k]==ppbase[n][l]!==ppbase[n][j] #ABAA
                spcount[6]=spcount[6]+counts[n]
            elseif ppbase[n][i]==ppbase[n][k]!==ppbase[n][j]==ppbase[n][l] #ABAB
                spcount[7]=spcount[7]+counts[n]
            elseif ppbase[n][i]==ppbase[n][k]!==ppbase[n][j]!==ppbase[n][l] #ABAC
                spcount[8]=spcount[8]+counts[n]
            elseif ppbase[n][i]==ppbase[n][l]!==ppbase[n][j]==ppbase[n][k] #ABBA
                spcount[9]=spcount[9]+counts[n]
            elseif ppbase[n][i]!==ppbase[n][j]==ppbase[n][k]==ppbase[n][l] #BAAA
                spcount[10]=spcount[10]+counts[n]
            elseif ppbase[n][i]!==ppbase[n][j]==ppbase[n][k]!==ppbase[n][l]!==ppbase[n][i] #ABBC
                spcount[11]=spcount[11]+counts[n]
            elseif ppbase[n][i]==ppbase[n][l]!==ppbase[n][j]!==ppbase[n][k] #CABC
                spcount[12]=spcount[12]+counts[n]
            elseif ppbase[n][j]==ppbase[n][l]!==ppbase[n][i]!==ppbase[n][k] #BACA
                spcount[13]=spcount[13]+counts[n]
            elseif ppbase[n][i]!==ppbase[n][j]!==ppbase[n][k]==ppbase[n][l] #BCAA
                spcount[14]=spcount[14]+counts[n]
            elseif ppbase[n][i]!==ppbase[n][j]!==ppbase[n][k]!==ppbase[n][l] #ABCD
                spcount[15]=spcount[15]+counts[n]
            else
                print("Something went wrong while getting the site pattern frequences... cannot find corresponding site pattern at site $(n)")
            end
        end
        push!(p.spcounts,spcount)
    end 
end

function spRearrange(p::Phylip)
    Allquartet=p.allquartet
    numQuarts=length(Allquartet)
    SPcounts=p.spcounts
    for n in 1:numQuarts
        i=Allquartet[n][1]
        j=Allquartet[n][2]
        k=Allquartet[n][3]
        l=Allquartet[n][4]

        AAAA=SPcounts[n][1]
        AAAB=SPcounts[n][2]
        AABA=SPcounts[n][3]
        AABB=SPcounts[n][4]
        AABC=SPcounts[n][5]
        ABAA=SPcounts[n][6]
        ABAB=SPcounts[n][7]
        ABAC=SPcounts[n][8]
        ABBA=SPcounts[n][9]
        BAAA=SPcounts[n][10]
        ABBC=SPcounts[n][11]
        CABC=SPcounts[n][12]
        BACA=SPcounts[n][13]
        BCAA=SPcounts[n][14]
        ABCD=SPcounts[n][15]
    
        #rearrage ijkl r=[i,j,k,l]#push!(p.allquartet,r)#1    
        r=[i,j,l,k]; push!(p.allquartet,r)#2
        r=[i,k,j,l]; push!(p.allquartet,r)#3
        r=[i,k,l,j]; push!(p.allquartet,r)#4
        r=[i,l,j,k]; push!(p.allquartet,r)#5
        r=[i,l,k,j]; push!(p.allquartet,r)#6
        r=[j,i,k,l]; push!(p.allquartet,r)#7
        r=[j,i,l,k]; push!(p.allquartet,r)#8
        r=[j,k,i,l]; push!(p.allquartet,r)#9
        r=[j,k,l,i]; push!(p.allquartet,r)#10
        r=[j,l,i,k]; push!(p.allquartet,r)#11
        r=[j,l,k,i]; push!(p.allquartet,r)#12
        r=[k,i,j,l]; push!(p.allquartet,r)#13
        r=[k,i,l,j]; push!(p.allquartet,r)#14
        r=[k,j,i,l]; push!(p.allquartet,r)#15
        r=[k,j,l,i]; push!(p.allquartet,r)#16
        r=[k,l,i,j]; push!(p.allquartet,r)#17
        r=[k,l,j,i]; push!(p.allquartet,r)#18
        r=[l,i,j,k]; push!(p.allquartet,r)#19
        r=[l,i,k,j]; push!(p.allquartet,r)#20
        r=[l,j,i,k]; push!(p.allquartet,r)#21
        r=[l,j,k,i]; push!(p.allquartet,r)#22
        r=[l,k,i,j]; push!(p.allquartet,r)#23
        r=[l,k,j,i]; push!(p.allquartet,r)#24

        #rearrage spc=[AAAA,AAAB,thrre...ABCD]#push!(p.spcounts,spc)#1
        spc=[AAAA,AABA,AAAB,AABB,AABC,ABAA,ABBA,CABC,ABAB,BAAA,BACA,ABAC,ABBC,BCAA,ABCD]; push!(p.spcounts,spc)#2
        spc=[AAAA,AAAB,ABAA,ABAB,ABAC,AABA,AABB,AABC,ABBA,BAAA,ABBC,CABC,BCAA,BACA,ABCD]; push!(p.spcounts,spc)#3
        spc=[AAAA,ABAA,AAAB,ABAB,ABAC,AABA,ABBA,CABC,AABB,BAAA,BCAA,AABC,ABBC,BACA,ABCD]; push!(p.spcounts,spc)#4
        spc=[AAAA,AABA,ABAA,ABBA,CABC,AAAB,AABB,AABC,ABAB,BAAA,BACA,ABAC,BCAA,ABBC,ABCD]; push!(p.spcounts,spc)#5
        spc=[AAAA,ABAA,AABA,ABBA,CABC,AAAB,ABAB,ABAC,AABB,BAAA,BCAA,AABC,BACA,ABBC,ABCD]; push!(p.spcounts,spc)#6
        spc=[AAAA,AAAB,AABA,AABB,AABC,BAAA,ABBA,ABBC,ABAB,ABAA,ABAC,BACA,CABC,BCAA,ABCD]; push!(p.spcounts,spc)#7    
        spc=[AAAA,AABA,AAAB,AABB,AABC,BAAA,ABAB,BACA,ABBA,ABAA,CABC,ABBC,ABAC,BCAA,ABCD]; push!(p.spcounts,spc)#8
        spc=[AAAA,AAAB,BAAA,ABBA,ABBC,AABA,AABB,AABC,ABAB,ABAA,ABAC,BACA,BCAA,CABC,ABCD]; push!(p.spcounts,spc)#9   
        spc=[AAAA,BAAA,AAAB,ABBA,ABBC,AABA,ABAB,BACA,AABB,ABAA,BCAA,AABC,ABAC,CABC,ABCD]; push!(p.spcounts,spc)#10
        spc=[AAAA,AABA,BAAA,ABAB,BACA,AAAB,AABB,AABC,ABBA,ABAA,CABC,ABBC,BCAA,ABAC,ABCD]; push!(p.spcounts,spc)#11
        spc=[AAAA,BAAA,AABA,ABAB,BACA,AAAB,ABBA,ABBC,AABB,ABAA,BCAA,AABC,CABC,ABAC,ABCD]; push!(p.spcounts,spc)#12
        spc=[AAAA,AAAB,ABAA,ABAB,ABAC,BAAA,ABBA,ABBC,AABB,AABA,AABC,BCAA,CABC,BACA,ABCD]; push!(p.spcounts,spc)#13
        spc=[AAAA,ABAA,AAAB,ABAB,ABAC,BAAA,AABB,BCAA,ABBA,AABA,CABC,ABBC,AABC,BACA,ABCD]; push!(p.spcounts,spc)#14
        spc=[AAAA,AAAB,BAAA,ABBA,ABBC,ABAA,ABAB,ABAC,AABB,AABA,AABC,BCAA,BACA,CABC,ABCD]; push!(p.spcounts,spc)#15
        spc=[AAAA,BAAA,AAAB,ABBA,ABBC,ABAA,AABB,BCAA,ABAB,AABA,BACA,ABAC,AABC,CABC,ABCD]; push!(p.spcounts,spc)#16
        spc=[AAAA,ABAA,BAAA,AABB,BCAA,AAAB,ABAB,ABAC,ABBA,AABA,CABC,ABBC,BACA,AABC,ABCD]; push!(p.spcounts,spc)#17
        spc=[AAAA,BAAA,ABAA,AABB,BCAA,AAAB,ABBA,ABBC,ABAB,AABA,BACA,ABAC,CABC,AABC,ABCD]; push!(p.spcounts,spc)#18
        spc=[AAAA,AABA,ABAA,ABBA,CABC,BAAA,ABAB,BACA,AABB,AAAB,AABC,BCAA,ABAC,ABBC,ABCD]; push!(p.spcounts,spc)#19
        spc=[AAAA,ABAA,AABA,ABBA,CABC,BAAA,AABB,BCAA,ABAB,AAAB,ABAC,BACA,AABC,ABBC,ABCD]; push!(p.spcounts,spc)#20
        spc=[AAAA,AABA,BAAA,ABAB,BACA,ABAA,ABBA,CABC,AABB,AAAB,AABC,BCAA,ABBC,ABAC,ABCD]; push!(p.spcounts,spc)#21     
        spc=[AAAA,BAAA,AABA,ABAB,BACA,ABAA,AABB,BCAA,ABBA,AAAB,ABBC,CABC,AABC,ABAC,ABCD]; push!(p.spcounts,spc)#22
        spc=[AAAA,ABAA,BAAA,AABB,BCAA,AABA,ABBA,CABC,ABAB,AAAB,ABAC,BACA,ABBC,AABC,ABCD]; push!(p.spcounts,spc)#23
        spc=[AAAA,BAAA,ABAA,AABB,BCAA,AABA,ABAB,BACA,ABBA,AAAB,ABBC,CABC,ABAC,AABC,ABCD]; push!(p.spcounts,spc)#24
    end
end

function spRearrangeCHECK(inputfile::AbstractString,display=false::Bool)
    p=Phylip(inputfile)
    UniqueBase,BaseCounts=PhylipFileInfo(inputfile, p)
    getUniqueQuartets(p)
    sitePatternCounts(p,UniqueBase,BaseCounts)
    laq0=length(p.allquartet)
    lsp0=length(p.spcounts)
    spRearrange(p)
    laq1=length(p.allquartet)
    lsp1=length(p.spcounts)
    println("
There were $laq0 quartets and $lsp0 site pattern frequenciesy before the function spRearrange,
and now there are $laq1 quartets and $lsp1 site pattern frequencies after excuting the function spRearrange.
The two numbers (number of quartes and site pattern frequencies)
Set [display=true] to see both quartets and coresponding site pattern frequencies results.")
    if(display)
        println("p.allquartet: $(p.allquartet)")
        println("p.spcounts: $(p.spcounts)")
    end
end

"""
    writeSitePatternCounts(p::Phylip,write::Bool,inputfile::AbstractString)
    writeSitePatternCounts(p::Phylip,inputfile::AbstractString)
    writeSitePatternCounts(inputfile::AbstractString)
    writeSitePatternCounts(p::Phylip)
    

Writes the observed site pattern frequencies from the phylip file as .csv file in the working 
directory and name it as sitePatternCounts_inputfile.csv. 
Although it requires an attribute `inputfile`, it does not actually read the file as long as
the Phylip object in memory is provided. (i.e., any `AbstractString` can be provided here that 
is used as a suffix of the ) However, if only the input file name is provided, 
it runs `readPhylipFile!` first then write csv. file.

### Input
`p` A Type Phylip object\\
`write`  Boolean variable to export site pattern frequences in .csv\\
`inputfile` Name of phylip file as a AbstractString\\
"""
function writeSitePatternCounts(p::Phylip,write::Bool,csvname::AbstractString,inputfile::AbstractString)
    
    df=DataFrame(i=Any[],j=Any[],k=Any[],l=Any[], 
                AAAA=Int[], AAAB=Int[], AABA=Int[], AABB=Int[], AABC=Int[], 
                ABAA=Int[], ABAB=Int[], ABAC=Int[], ABBA=Int[], BAAA=Int[], 
                ABBC=Int[], CABC=Int[], BACA=Int[], BCAA=Int[], ABCD=Int[])
    
    for n in 1:length(p.allquartet)
        firsttaxa   = p.allquartet[n][1]
        secondtaxa  = p.allquartet[n][2]
        thirdtaxa   = p.allquartet[n][3]
        forthtaxa   = p.allquartet[n][4]
        
        push!(df.i,p.nametaxa[firsttaxa])
        push!(df.j,p.nametaxa[secondtaxa])
        push!(df.k,p.nametaxa[thirdtaxa])
        push!(df.l,p.nametaxa[forthtaxa])

        push!(df.AAAA,p.spcounts[n][1])
        push!(df.AAAB,p.spcounts[n][2])
        push!(df.AABA,p.spcounts[n][3])
        push!(df.AABB,p.spcounts[n][4])
        push!(df.AABC,p.spcounts[n][5])
        push!(df.ABAA,p.spcounts[n][6])
        push!(df.ABAB,p.spcounts[n][7])
        push!(df.ABAC,p.spcounts[n][8])
        push!(df.ABBA,p.spcounts[n][9])
        push!(df.BAAA,p.spcounts[n][10])
        push!(df.ABBC,p.spcounts[n][11])
        push!(df.CABC,p.spcounts[n][12])
        push!(df.BACA,p.spcounts[n][13])
        push!(df.BCAA,p.spcounts[n][14])
        push!(df.ABCD,p.spcounts[n][15])
    end
    if write==true 
        if isempty(csvname)
            CSV.write("sitePatternCounts_$inputfile.csv",df)
            println("A [.csv] file is saved as sitePatternCounts_$inputfile.csv in the current working directory.")
        else
            CSV.write("$csvname.csv",df)
            println("A [.csv] file is saved as $csvname.csv in the current working directory.")
        end
    end
end
#function writeSitePatternCounts(inputfile::AbstractString) readPhylipFile!(inputfile,true) end
writeSitePatternCounts(p::Phylip,csvname::AbstractString)=writeSitePatternCounts(p,true,csvname,csvname)
#writeSitePatternCounts(p::Phylip,inputfile::AbstractString)=writeSitePatternCounts(p,true,inputfile)

function binaryIndexforQuartet(p::Phylip)
    Allquartet=p.allquartet
    numQuarts=length(Allquartet)
    for n in 1:numQuarts
        i=bitstring(Int8(Allquartet[n][1]))
        j=bitstring(Int8(Allquartet[n][2]))
        k=bitstring(Int8(Allquartet[n][3]))
        l=bitstring(Int8(Allquartet[n][4]))
        push!(p.index,[i*j*k*l])
    end
end

function printCheckPoint(p::Phylip)
    open("$(p.filename).ckp", "w") do file
        write(file, "$(p.filename)\n")
        write(file, "$(p.time)\n")
        write(file, "$(p.numtaxa)\n")
        write(file, "$(p.seqleng)\n")
        
        for element in p.nametaxa
            write(file, "$(element),")
        end
        write(file, "\n")
        
        for element in p.counttaxa
            write(file, "$(element),")
        end
        write(file, "\n")
        
        for element in p.allquartet
            for leaves in element
                write(file, "$(leaves),")
            end
            write(file, ";")
        end
        write(file, "\n")
        
        for element in p.index
            for index in element
                write(file, "$(index),")    
            end
        end
        write(file, "\n")

        for spc in p.spcounts
            for eachpattern in spc
                write(file, "$(eachpattern),")
            end
            write(file, ";")
        end
    end
end


"""
    readCheckPoint(ckpfile::AbstractString)

Reads in `.ckp` file. gegeege



"""
function readCheckPoint(ckpfile::AbstractString)
    p=Phylip()
    nline=0
    for line in eachline(ckpfile)
        nline+=1
        if nline==1; p.filename=line; end
        if nline==2; p.time=parse(Float64,line); end
        if nline==3; p.numtaxa=parse(Float64,line); end
        if nline==4; p.seqleng=parse(Float64,line); end
        if nline==5; 
            line=chop(line,tail=1)
            nametaxa=split(line, ",")
            for element in nametaxa
                push!(p.nametaxa,element)
            end
        end
        if nline==6; 
            line=chop(line,tail=1)
            counttaxa=split(line, ",")         
            for element in counttaxa
                element=parse(Int64, element)
                push!(p.counttaxa,element)
            end
        end
        if nline==7; 
            line=chop(line,tail=1)
            allquarts=split(line, ";")     
            for eachquart in allquarts
                eachquart=chop(eachquart,tail=1)
                eachquart=split(eachquart, ",")
                leaf1=parse(Int64, eachquart[1])
                leaf2=parse(Int64, eachquart[2])
                leaf3=parse(Int64, eachquart[3])
                leaf4=parse(Int64, eachquart[4])
                push!(p.allquartet,[leaf1,leaf2,leaf3,leaf4])
            end
        end
        if nline==8;
            line=chop(line,tail=1)
            index=split(line, ",")
            for eachindex in index
                ind=[eachindex]
                push!(p.index,ind)
            end
        end
        if nline==9;
            line=chop(line,tail=1)
            spc=split(line,";")
            for eachpattern in spc
                eachpattern=chop(eachpattern,tail=1)
                eachpattern=split(eachpattern, ",")
                pattern1=parse(Float64, eachpattern[1])
                pattern2=parse(Float64, eachpattern[2])
                pattern3=parse(Float64, eachpattern[3])
                pattern4=parse(Float64, eachpattern[4])
                pattern5=parse(Float64, eachpattern[5])
                pattern6=parse(Float64, eachpattern[6])
                pattern7=parse(Float64, eachpattern[7])
                pattern8=parse(Float64, eachpattern[8])
                pattern9=parse(Float64, eachpattern[9])
                pattern10=parse(Float64, eachpattern[10])
                pattern11=parse(Float64, eachpattern[11])
                pattern12=parse(Float64, eachpattern[12])
                pattern13=parse(Float64, eachpattern[13])
                pattern14=parse(Float64, eachpattern[14])
                pattern15=parse(Float64, eachpattern[15])
                push!(p.spcounts,[pattern1,pattern2,pattern3,pattern4,pattern5,pattern6,pattern7,pattern8,pattern9,pattern10,pattern11,pattern12,pattern13,pattern14,pattern15])
            end
        end

    end
    return p
end


function readCSVFile(inputfile::AbstractString)
    df=DataFrame(CSV.File(inputfile))
    return df
end