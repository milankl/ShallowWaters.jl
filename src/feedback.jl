function readable_secs(secs::Real)
    #= Returns a human readable string representing seconds in terms of days, hours, minutes, seconds. =#

    days = Int(floor(secs/3600/24))
    hours = Int(floor((secs/3600) % 24))
    minutes = Int(floor((secs/60) % 60))
    seconds = Int(floor(secs%3600%60))

    if days > 0
        return "$(days)d, $(hours)h"
    elseif hours > 0
        return "$(hours)h, $(minutes)min"
    elseif minutes > 0
        return "$(minutes)min, $(seconds)s"
    else
        return "$(seconds)s"
    end
end

function duration_estimate(i::Int,t::Real,nt::Int)
    #= Estimates the total time the model integration will take.=#
    time_per_step = (time()-t) / (i-10)
    time_total = Int(round(time_per_step*nt))
    time_to_go = Int(round(time_per_step*(nt-i)))
    println("Model integration will take approximately "*readable_secs(time_total)*",")
    println("and is hopefully done on "*Dates.format(now() + Dates.Second(time_to_go),Dates.RFC1123Format))
end

function feedback_ini()
    println("Starting juls on "*Dates.format(now(),Dates.RFC1123Format))
    return time()
end

function feedback_end(t::Real)
    println("Integration done in "*readable_secs(time()-t)*".")
end

function feedback(i::Int,t::Real,nt::Int)
    if i == 10
        return time()    # measure time after 10 loops to avoid overhead
    elseif i == 100
        duration_estimate(i,t,nt)
        return t
    else
        return t
    end
end
