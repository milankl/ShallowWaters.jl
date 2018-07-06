const ω = 2π/(24.*3600.)    # Earth's angular frequency [s^-1]
const R = 6.371e6           # Earth's radius [m]
const f_0 = 2*ω*sind(ϕ)     # Coriolis parameter [s^-1]
const β = 2*ω/R*cosd(ϕ)     # Coriolis derivative wrt latitude [(ms)^-1]

function beta_plane()

    xx_u,yy_u = meshgrid(x_u,y_u)
    xx_v,yy_v = meshgrid(x_v,y_v)
    xx_q,yy_q = meshgrid(x_q,y_q)

    f_u = Numtype.(Δ*(f_0 + β*yy_u))
    f_v = Numtype.(Δ*(f_0 + β*yy_v))
    f_q = Numtype.(Δ*(f_0 + β*yy_q))
    #f_u = Numtype.(f_0 + β*yy_u)
    #f_v = Numtype.(f_0 + β*yy_v)
    #f_q = Numtype.(f_0 + β*yy_q)
    return f_u,f_v,f_q
end
