function greet()
    return Dates.now()
end

function checkDEBUG(i::Any)
    println(i)
    @debug "Check debug logging"
end
