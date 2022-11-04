function greet()
    now=Dates.now()
    msg=println("
    Thanks for using PhyNEST!
    It seems like PhyNE has been loaded correctly as of $now.
    Let's make some networks...:)! 
    [Sungsik Kong, October 2022, @MBD OSU]")
    return msg
end

function checkDEBUG(i::Any)
    println(i)
    @debug "Check debug logging"
end
