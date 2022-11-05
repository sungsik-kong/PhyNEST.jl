function greet()
    now=Dates.now()
    msg=println("
    Thank you for using PhyNEST!
    Current time is: $now.")
    return msg
end

function checkDEBUG(i::Any)
    println(i)
    @debug "Check debug logging"
end
