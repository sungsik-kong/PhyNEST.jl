function greet()
    now=Dates.now()
    msg=println("Thank you for using PhyNEST: Phylogenetic Network Estimation using SiTe patterns.
Please report bugs or make suggestions to https://groups.google.com/g/phynest-users.
Current time is $now.")
    return msg
end

function checkDEBUG(i::Any)
    println(i)
    @debug "Check debug logging"
end
