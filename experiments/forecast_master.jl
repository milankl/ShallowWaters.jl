using Juls

RunJuls(Float32,
        output=true,
        Ndays=500.0,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast",
        output_vars=["u","v","η"],
        α=1.,
        diffusion="constant")

RunJuls(Float32,
        output=true,
        Ndays=50*365,
        outpath="/network/aopp/chaos/pred/kloewer/julsdata/forecast",
        output_vars=["u","v","η"],
        initial_cond="ncfile",
        init_run_id=0,
        α=1.,
        diffusion="constant")
