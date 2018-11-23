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

function duration_estimate(i,t,nt,progrtxt)
    #= Estimates the total time the model integration will take.=#
    time_per_step = (time()-t) / (i-50)
    time_total = Int(round(time_per_step*nt))
    time_to_go = Int(round(time_per_step*(nt-i)))

    s1 = "Model integration will take approximately "*readable_secs(time_total)*","
    s2 = "and is hopefully done on "*Dates.format(now() + Dates.Second(time_to_go),Dates.RFC1123Format)

    println(s1)     # print inline
    println(s2)
    if output == 1  # print in txt
        write(progrtxt,"\n"*s1*"\n")
        write(progrtxt,s2*"\n")
        flush(progrtxt)
    end
end

function nan_detection(u::AbstractMatrix,v::AbstractMatrix,η::AbstractMatrix)
    #= Returns a boolean  input matrices u,v,η contains a NaN.
    TODO include a check for Posits
    =#
    n_nan = sum(isnan.(u)) + sum(isnan.(v)) + sum(isnan.(η))
    if n_nan > 0
        return true
    else
        return false
    end

    #TODO include of check for tracer!

end

function progress_txt_ini()
    # initialises a txt file for progress feedback
    progrtxt = open(runpath*"progress.txt","w")
    return progrtxt
end

function feedback_ini()
    # initialising the progress txt file in case of output
    if output == 1
        progrtxt = progress_txt_ini()
        s = "Starting juls run $run_id on "*Dates.format(now(),Dates.RFC1123Format)
        println(s)
        write(progrtxt,s*"\n")
        write(progrtxt,"Juls will integrate $(Ndays)days at a resolution of $(nx)x$(ny) with Δ=$(Δ/1e3)km\n")
        write(progrtxt,"Initial conditions are ")
        if initial_cond == "rest"
            write(progrtxt,"rest.\n")
        else
            write(progrtxt,"last time step of run $init_run_id.\n")
        end
        write(progrtxt,"Boundary conditions are $bc_x with lbc=$lbc.\n")
        write(progrtxt,"Numtype is "*string(Numtype)*".\n")
        write(progrtxt,"\nAll data will be stored in $runpath\n")
    else
        println("Starting juls on "*Dates.format(now(),Dates.RFC1123Format))
        progrtxt = nothing
    end

    return time(),progrtxt
end

function feedback(u,v,η,sst,i,t,nt,nans_detected,progrtxt)
    if i == 50
        t = time()    # measure time after 50 loops to avoid overhead and make sure tracer advection executed once
    elseif i == 100
        duration_estimate(i,t,nt,progrtxt)
    end

    if !nans_detected
        if i % nout == 0    # only check for nans when output is produced
            nans_detected = nan_detection(u,v,η)
            if nans_detected
                println(" NaNs detected at time step $i")
                if output == 1
                    write(progrtxt," NaNs detected at time step $i")
                    flush(progrtxt)
                end
            end
        end
    end

    if i > 100      # show percentage only after duration is estimated
        progress(i,nt,progrtxt)
    end

    return t,nans_detected
end

function feedback_end(progrtxt,t::Real)
    # finalise the progress txt file.
    s = " Integration done in "*readable_secs(time()-t)*"."
    println(s)
    if output == 1
        write(progrtxt,"\n"*s[2:end]*"\n")  # close txt file with last output
        flush(progrtxt)
    end
end

function progress(i,nt,progrtxt)
    # progress feedback in percent
    if ((i+1)/nt*100 % 1) < (i/nt*100 % 1)  # update every 1 percent steps.
        percent = Int(round((i+1)/nt*100))
        print("\r\u1b[K")
        print("$percent%")
        if (output == 1) && (percent % 5 == 0)
            write(progrtxt,"\n$percent%")
            flush(progrtxt)
        end
    end
end
