include("../Plots.jl")
include("../Check_Monotonicity.jl")
include("../Optimizations.jl")
include("../Interpolations.jl")
include("../Parameters.jl")
include("../1_Income.jl")
include("../2_Learning.jl")
include("NHL_3_Grid.jl")
include("NHL_4_Retirement.jl")
include("NHL_5_Transition.jl")
include("NHL_6_Bellman.jl")
include("NHL_7_Simulate.jl")

# 1. Draw Income Distribution
guvenen_distribution = true

if guvenen_distribution
  (yit, α, β, ymedian, pension) =
    incomeDistribution(agents, bs,
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")
else
  (yit, α, β, ymedian, pension) =
    incomeDistribution(agents, bs, μₐ, μᵦ, var_a, var_b, var_ɛ, ρ, var_η,
                       br, tW)
end

# 2. Construct individual specific belief histories
(s_f_i, stdy) = learning(agents, bs, tW, α, β, yit, ρ,
                         var_η, var_ɛ, guvenen_distribution)
# 3. Construct Grids
(wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, stdy, wpoints, apoints, bpoints, zpoints,
        wpoints_R, ypoints_R, wmaxR, power, r, tR, false, true)

# 4. Solve Retirement Problem
(v_R, wp_R) =
  solveRetirement(wgrid_R, ygrid_R, wpoints_R, ypoints_R, r, δ)

# 5. Solve Transition Problem
#(v, wp, constrained) =
#  solveTransition(v_R, wgrid_R, ygrid_R, wgrid, agrid, bgrid, ymedian, r, δ, tW)

# 5.1 Solve Terminal Period Problem instead of 4. and 5.
(wp, c, v, c_over_x) =
  solveTransition(wgrid, agrid, bgrid, zgrid, r, δ)

# 6. Solve Working Life Problem
(v, wp, c_over_x) =
  solveWorkingLife(v, wp, c_over_x, wgrid, agrid, bgrid, zgrid, stdy, r, δ)

# 7. Simulate consumption and wealth distributions
(c_t, w_t, wp_t) =
  simulate(wp, wgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)
