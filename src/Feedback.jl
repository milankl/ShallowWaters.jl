"""Returns a human readable string representing seconds in terms of days, hours, minutes or seconds."""
function readable_secs(secs::Real)
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

"""Estimates the total time the model integration will take."""
function duration_estimate(i,t,nt,progrtxt)
    time_per_step = (time()-t) / (i-nadvstep)
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

"""Returns a boolean whether the prognostic variables contains a NaN."""
function nan_detection(u::AbstractMatrix,v::AbstractMatrix,η::AbstractMatrix,sst::AbstractMatrix)
    #TODO include a check for Posits, are posits <: AbstractFloat?
    #TODO include check for tracer by other means than nan? (semi-Lagrange is unconditionally stable...)

    n_nan = sum(isnan.(u)) + sum(isnan.(v)) + sum(isnan.(η)) + sum(isnan.(sst))
    if n_nan > 0
        return true
    else
        return false
    end
end

"""Initialises the progress txt file."""
function feedback_ini()
    if output
        progrtxt = open(runpath*"progress.txt","w")
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
        write(progrtxt,"Time steps are (Lin,Diff,Advcor,Lagr,Output)\n")
        write(progrtxt,"$dtint, $(dtint*nstep_diff), $(dtint*nstep_advcor), $dtadvint, $(output_dt*3600)\n")
        write(progrtxt,"\nAll data will be stored in $runpath\n")
    else
        println("Starting juls on "*Dates.format(now(),Dates.RFC1123Format)*" without output.")
        progrtxt = nothing
    end

    return time(),progrtxt
end

"""Feedback function that calls duration estimate, nan_detection and progress."""
function feedback(u,v,η,sst,i,t,nt,nans_detected,progrtxt)
    if i == nadvstep # measure time after tracer advection executed once
        t = time()
    elseif i == min(2*nadvstep,nadvstep+50)
        # after the tracer advection executed twice or at least 50 steps
        duration_estimate(i,t,nt,progrtxt)
    end

    if !nans_detected
        if i % nout == 0    # only check for nans when output is produced
            nans_detected = nan_detection(u,v,η,sst)
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

"""Finalises the progress txt file."""
function feedback_end(progrtxt,t::Real)
    s = " Integration done in "*readable_secs(time()-t)*"."
    println(s)
    if output
        write(progrtxt,"\n"*s[2:end]*"\n")  # close txt file with last output
        flush(progrtxt)
    end
end

"""Converts time step into percent for feedback."""
function progress(i,nt,progrtxt)
    if ((i+1)/nt*100 % 1) < (i/nt*100 % 1)  # update every 1 percent steps.
        percent = Int(round((i+1)/nt*100))
        print("\r\u1b[K")
        print("$percent%")
        if output && (percent % 5 == 0) # write out only every 5 percent step.
            write(progrtxt,"\n$percent%")
            flush(progrtxt)
        end
    end
end
