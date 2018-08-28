# Earth's constants
const ω = 2π/(24*3600)    # Earth's angular frequency [s^-1]
const R = 6.371e6           # Earth's radius [m]
const f₀ = 2*ω*sind(ϕ)      # Coriolis parameter [s^-1]
const β = 2*ω/R*cosd(ϕ)     # Coriolis derivative wrt latitude [(ms)^-1]

function beta_plane()

    if bc_x == "periodic"
        # points on the right edge needed too
        xx_q,yy_q = meshgrid(x_q_halo[3:end-1],y_q_halo[3:end-2])
    else
        xx_q,yy_q = meshgrid(x_q,y_q)
    end

    # for non-dimensional gradient operators f contains the grid spacing Δ
    f_q = Numtype.(Δ*(f₀ .+ β*yy_q))

    return f_q
end
