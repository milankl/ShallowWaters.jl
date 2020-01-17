include("/home/kloewer/git/Juls.jl/src/Juls.jl")
using .Juls
using FileIO

"""Finds the first gap in a list of integers."""
function gap(a::Array{Int,1})
    try
        return minimum([i for i in minimum(a):maximum(a) if ~(i in a)])
    catch
        return maximum(a)+1
    end
end

function get_run_id(path::String,order::String="continue")
    runlist = filter(x->startswith(x,"run"),readdir(path))
    existing_runs = [parse(Int,id[4:end]) for id in runlist]
    if length(existing_runs) == 0           # if no runfolder exists yet
        run_id = 0
    else                                    # create next folder
        if order == "fill"  # find the smallest gap in runfolders
            run_id = gap(existing_runs)
        elseif order == "continue" # find largest folder and count one up
            run_id = maximum(existing_runs)+1
        end
    end
    return run_id
end

path = "/network/aopp/chaos/pred/kloewer/julsdata/forecast/"
startis = load(joinpath(path,"starti.jld2"))["starti"]
outpath = joinpath(path,"Float64")

for i in 1:10
    run_id = get_run_id(outpath)

    starti = startis[run_id+1]

    RunJuls(Float64,
            output=true,
            Ndays=100.0,
            outpath=outpath,
            initial_cond="ncfile",
            output_vars=["u","v","η","sst"],
            initpath=path,
            init_run_id=1,
            init_starti=starti
            )

end
