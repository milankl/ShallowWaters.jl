# # ShallowWaters.jl - The equations
#
# The shallow water equations for the prognostic variables velocity $\mathbf{u} = (u,v)$ and sea surface elevation $\eta$ over the 2-dimensional domain $\Omega$ in $x,y$ are
#
# ```math
# \begin{align}
# \partial_t u &+ u\partial_xu + v\partial_yu - fv = -g\partial_x\eta + D_x(u,v,\eta) + F_x, \\
# \partial_t v &+ u\partial_xv + v\partial_yv + fu = -g\partial_y\eta + D_y(u,v,\eta) + F_y, \\
# \partial_t \eta &+ \partial_x(uh) + \partial_y(vh) = 0.
# \end{align}
# ```
# where the first two are the momentum equations for $u,v$ and the latter is the continuity equation. The layer thickness is $h = \eta + H$ with $H=H(x,y)$ being the bottom topography. The gravitational acceleration is $g$, the coriolis parameter $f = f(y)$ depends only (i.e. latitude) only and the beta-plane approximation is used. $(D_x,D_y) = \mathbf{D}$ are the dissipative terms
# 
# ```math
# \mathbf{D} = -\frac{c_D}{h}\vert \mathbf{u} \vert \mathbf{u} - u \nabla^4\mathbf{u}
# ```

# which are a sum of a quadratic bottom drag with dimensionless coefficient $c_D$ and a biharmonic diffusion with viscosity coefficient $u$. $\mathbf{F} = (F_x,F_y)$ is the wind forcing which can depend on time and space, i.e. $F_x = F_x(x,y,t)$ and $F_y = F_y(x,y,t)$.

# ## The vector invariant formulation

# The Bernoulli potential $p$ is introduced as

# ```math
# p = \frac{1}{2}(u^2 + v^2) + gh
# ```

# The relative vorticity $\zeta = \partial_xv + \partial_yu$ lets us define the potential vorticity $q$ as

# ```math
# q = \frac{f + \zeta}{h}
# ```

# such that we can rewrite the shallow water equations as

# ```math
# \begin{align}
# \partial_t u &= qhv -g\partial_xp + D_x(u,v,\eta) + F_x, \\
# \partial_t v &= -qhu -g\partial_yp + D_y(u,v,\eta) + F_y, \\ 
# \partial_t \eta &= -\partial_x(uh) -\partial_y(vh).
# \end{align}
# ```

# ## Runge-Kutta time discretisation

# Let
# ```math
# R(u,v,\eta) = \begin{pmatrix} 
#                 qhv-g\partial_xp+F_x 
#                 -qhu-g\partial_yp+F_y 
#                 -\partial_x(uh) -\partial_y(vh)
#                 \end{pmatrix}
# ```
# 
# be the non-dissipative right-hand side, i.e. excluding the dissipative terms $\mathbf{D}$. Then we dicretise the time derivative with 4th order Runge-Kutta with $\mathbf{k}_n = (u_n,v_n,\eta_n)$ by
# 
# ```math
# \begin{align}
# \mathbf{d}_1 &= R(\mathbf{k}_n) , \\
# \mathbf{d}_2 &= R(\mathbf{k}_n + \frac{1}{2}\Delta t \mathbf{d}_1), \\
# \mathbf{d}_3 &= R(\mathbf{k}_n + \frac{1}{2}\Delta t \mathbf{d}_2), \\
# \mathbf{d}_4 &= R(\mathbf{k}_n + \Delta t \mathbf{d}_3), \\
# u_{n+1}^*,v_{n+1}^*,\eta_{n+1} = \mathbf{k}_{n+1} &= \mathbf{k}_n + \frac{1}{6}\Delta t(\mathbf{d}_1 + 2\mathbf{d}_2 + 2\mathbf{d}_3 + \mathbf{d}_1),
# \end{align}
# ```

# and the dissipative terms are then added semi-implictly.

# ```math
# \begin{align}
# u_{n+1} = u_{n+1}^* + \Delta t D_x(u_{n+1}^*,v_{n+1}^*,\eta_{n+1}) , \\
# v_{n+1} = v_{n+1}^* + \Delta t D_y(u_{n+1}^*,v_{n+1}^*,\eta_{n+1}).
# \end{align}
# ```


# Consequently, the dissipative terms only have to be evaluated once per time step, which reduces the computational cost of the right=hand side drastically. This is motivated as the Courant number $C = \sqrt{gH_0}$ is restricted by gravity waves, which are caused by $\partial_t\mathbf{u} = -g\nabla\eta$ and $\partial_t\eta = -H_0\nabla \cdot \mathbf{u}$. The other terms have longer time scales, and it is therefore sufficient to solve those with larger time steps. With this scheme, the shallow water model runs stable at $C=1$."

# ## Strong stability preserving Runge-Kutta with semi-implicit continuity equation

# We split the right-hand side into the momentum equations and the continuity equation

# ```math
# R_m(u,v,\eta) = \begin{pmatrix} 
#                 qhv-g\partial_xp+F_x 
#                 -qhu-g\partial_yp+F_y 
#                 \end{pmatrix}, \quad R_\eta(u,v,\eta) = -\partial_x(uh) -\partial_y(vh)
# ```

# The 4-stage strong stability preserving Runge-Kutta scheme, with a semi-implicit treatment of the continuity equation then reads as

# ```math
# \begin{align}
# \mathbf{u}_1 &= \mathbf{u}_n + \frac{1}{2} \Delta t R_m(u_n,v_n,\eta_n), \quad \text{then}
# \quad \eta_1 = \eta_n + \frac{1}{2} \Delta t R_\eta(u_1,v_1,\eta_n), \\
# \mathbf{u}_2 &= \mathbf{u}_1 + \frac{1}{2} \Delta t R_m(u_1,v_1,\eta_1), \quad \text{then}
# \quad \eta_2 = \eta_1 + \frac{1}{2} \Delta t R_\eta(u_2,v_2,\eta_1) , \\
# \mathbf{u}_3 &= \frac{2}{3}\mathbf{u}_n + \frac{1}{3}\mathbf{u}_2 + \frac{1}{6} \Delta t R_m(u_2,v_2,\eta_2), \quad \text{then}
# \quad \eta_3 = \frac{2}{3}\eta_n + \frac{1}{3}\eta_2 + \frac{1}{6} \Delta t R_\eta(u_3,v_3,\eta_2), \\
# \mathbf{u}_{n+1} &= \mathbf{u}_3 + \frac{1}{2} \Delta t R_m(u_3,v_3,\eta_3), \quad \text{then}
# \quad \eta_{n+1} = \eta_3 + \frac{1}{2} \Delta t R_\eta(u_{n+1},v_{n+1},\eta_3).
# \end{align}
# ```

# ## Splitting the continuity equation

# From
# ```math
# \partial_t = -\partial_x(uh) - \partial_y(vh)
# ```
