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

function nan_detection(u::Array,v::Array,η::Array)
    # TODO include a check for Posit, Integers?
    n_nan = sum(isnan.(u)) + sum(isnan.(v)) + sum(isnan.(η))
    if n_nan > 0
        return true
    else
        return false
    end
end

function feedback_ini()
    if output == 1
        println("Starting juls run $run_id on "*Dates.format(now(),Dates.RFC1123Format))
    else
        println("Starting juls on "*Dates.format(now(),Dates.RFC1123Format))
    end

    return time()
end

function feedback_end(t::Real)
    println(" Integration done in "*readable_secs(time()-t)*".")
end

function feedback(u::Array,v::Array,η::Array,i::Int,t::Real,nt::Int,nans_detected::Bool)
    if i == 10
        t = time()    # measure time after 10 loops to avoid overhead
    elseif i == 100
        duration_estimate(i,t,nt)
    end

    if !nans_detected
        if i % nout == 0    # only check for nans when output is produced
            nans_detected = nan_detection(u,v,η)
            if nans_detected
                println(" NaNs detected at time step $i")
            end
        end
    end

    if i > 100      # show percentage only after duration is estimated
        progress(i,nt)
    end

    return t,nans_detected
end

function progress(i::Int,nt::Int)
    if ((i+1)/nt*100 % 1) < (i/nt*100 % 1)  # update every 1 percent steps.
        percent = Int(round((i+1)/nt*100))
        print("\r\u1b[K")
        print("$percent%")
    end
end
