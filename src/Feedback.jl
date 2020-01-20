@with_kw mutable struct Feedback
    t0::Float64=time()                      # start time
    t1::Float64=time()                      # time for duration estimation
    nans_detected::Bool=false               # did NaNs occur in the simulation?
    progress_txt::Union{IOStream,Nothing}   # txt is a Nothing in case of no output
    output::Bool                            # output to netCDF?
    i::Int=0                                # time step increment
    nt::Int                                 # number of time steps
    nout::Int                               # number of time steps with output
    run_id::Int                             # run identification number
    runpath::String                         # output path plus run????/
end


"""Returns a human readable string representing seconds in terms of days, hours, minutes or seconds."""
function readable_secs(secs::Real)
    days = Int(floor(secs/3600/24))
    hours = Int(floor((secs/3600) % 24))
    minutes = Int(floor((secs/60) % 60))
    seconds = Int(floor(secs%3600%60))
    secs1f = @sprintf "%.1fs" secs%3600%60
    secs2f = @sprintf "%.2fs" secs%3600%60

    if days > 0
        return "$(days)d, $(hours)h"
    elseif hours > 0
        return "$(hours)h, $(minutes)min"
    elseif minutes > 0
        return "$(minutes)min, $(seconds)s"
    elseif seconds > 10
        return secs1f
    else
        return secs2f
    end
end

"""Estimates the total time the model integration will take."""
function duration_estimate(feedback::Feedback,S::ModelSetup)

    @unpack t1,i,nt,output,progress_txt = feedback
    @unpack nadvstep = S.grid

    time_per_step = (time()-t1) / (i-nadvstep)
    time_total = Int(round(time_per_step*nt))
    time_to_go = Int(round(time_per_step*(nt-i)))

    s1 = "Model integration will take approximately "*readable_secs(time_total)*","
    s2 = "and is hopefully done on "*Dates.format(now() + Dates.Second(time_to_go),Dates.RFC1123Format)

    print("\r\u1b[K")
    println(s1)     # print inline
    println(s2)
    if output == 1  # print in txt
        write(progress_txt,"\n"*s1*"\n")
        write(progress_txt,s2*"\n")
        flush(progress_txt)
    end
end

"""Returns a boolean whether the prognostic variables contains a NaN."""
function nan_detection!(Prog::PrognosticVars,feedback::Feedback)

    #TODO include a check for Posits
    #TODO include check for sst

    @unpack u,v,η,sst = Prog

    n_nan = sum(isnan.(u)) + sum(isnan.(v)) + sum(isnan.(η)) + sum(isnan.(sst))
    if n_nan > 0
        feedback.nans_detected = true
    end
end

"""Initialises the progress txt file."""
function feedback_init(S::ModelSetup)
    @unpack output,init_starti = S.parameters
    @unpack nt,nout = S.grid

    if output

        @unpack Ndays,initial_cond,init_run_id,T,bc,α = S.parameters
        @unpack nx,ny,Δ = S.grid

        run_id,runpath = get_run_id_path(S)

        txt = open(joinpath(runpath,"progress.txt"),"w")
        s = "Starting Juls run $run_id on "*Dates.format(now(),Dates.RFC1123Format)
        println(s)
        write(txt,s*"\n")
        write(txt,"Juls will integrate $(Ndays)days at a resolution of $(nx)x$(ny) with Δ=$(Δ/1e3)km\n")
        write(txt,"Initial conditions are ")
        if initial_cond == "rest"
            write(txt,"rest.\n")
        else
            if init_starti > 0
                write(txt,"time step $init_starti of run $init_run_id.\n")
            else
                write(txt,"last time step of run $init_run_id.\n")
            end
        end
        write(txt,"Boundary conditions are $bc with lbc=$α.\n")
        write(txt,"Number format is "*string(T)*".\n")
        write(txt,"\nAll data will be stored in $runpath\n")

        # Parameter.txt
        ptxt = open(joinpath(runpath,"parameter.txt"),"w")
        print(ptxt,S.parameters)
        close(ptxt)
    else
        println("Starting Juls on "*Dates.format(now(),Dates.RFC1123Format)*" without output.")
        txt = nothing
        run_id = -1
        runpath = ""
    end

    return Feedback(progress_txt=txt,output=output,nt=nt,nout=nout,run_id=run_id,runpath=runpath)
end

"""Feedback function that calls duration estimate, nan_detection and progress."""
function feedback!(Prog::PrognosticVars,feedback::Feedback,S::ModelSetup)

    @unpack nadvstep = S.grid
    @unpack i = feedback

    if feedback.output
        if i == nadvstep # measure time after tracer advection executed once
            feedback.t1 = time()
        elseif i == 2*nadvstep
            # after the tracer advection executed twice or at least 50 steps
            duration_estimate(feedback,S)
        end
    end

    @unpack i, nans_detected, nout, output = feedback

    if !nans_detected
        if i % nout == 0    # only check for nans when output is produced
            nan_detection!(Prog,feedback)
            if feedback.nans_detected
                println(" NaNs detected at time step $i")
                if output
                    write(progress_txt," NaNs detected at time step $i")
                    flush(progress_txt)
                end
            end
        end
    end

    if feedback.i > 100      # show percentage only after duration is estimated
        progress!(feedback)
    end
end

"""Finalises the progress txt file."""
function feedback_end!(feedback::Feedback)
    @unpack output,t0,progress_txt = feedback

    s = " Integration done in "*readable_secs(time()-t0)*"."
    println(s)
    if output
        write(progress_txt,"\n"*s[2:end]*"\n")  # close txt file with last output
        flush(progress_txt)
    end
end

"""Converts time step into percent for feedback."""
function progress!(feedback::Feedback)

    @unpack i,nt,progress_txt,output = feedback

    if ((i+1)/nt*100 % 1) < (i/nt*100 % 1)  # update every 1 percent steps.
        percent = Int(round((i+1)/nt*100))
        print("\r\u1b[K")
        print("$percent%")
        if output && (percent % 5 == 0) # write out only every 5 percent step.
            write(progress_txt,"\n$percent%")
            flush(progress_txt)
        end
    end
end
