include("../Plots.jl")
include("../Chk_Monot.jl")
include("../Optimizations.jl")
include("../Interpolations.jl")
include("../Parameters.jl")
include("../1_Income.jl")
include("../2_Learning.jl")
include("HL_3_Grid.jl")
include("HL_4_Retirement.jl")
include("HL_5_Transition.jl")
include("HL_6_Bellman.jl")
include("HL_7_Simulate.jl")

# 1. Draw Income Distribution
guvenen_distribution = true

if guvenen_distribution
  (yit, α, β, ymedian, pension) =
    incomeDistribution(
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")
else
  (yit, α, β, ymedian, pension) =
    incomeDistribution(agents, bs, μₐ, μᵦ, var_a, var_b, var_ɛ, ρ, var_η,
                       br, tW)
end

# 2. Construct individual specific belief histories
(s_f_i, stdy) = learning(α, β, yit, ρ,
                         var_η, var_ɛ, guvenen_distribution)
# 3. Construct Grids
(wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R) =
    grids(s_f_i, stdy, wpoints, hpoints, apoints, bpoints, zpoints,
          wpoints_R, hpoints_R, ypoints_R, wmaxR, power, r, tR, false, true)

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, hgrid_R, ygrid_R, r, δ, λ)

# 5. Solve Transition Problem
(v, wp) = solveTransition(v_R, wgrid_R, hgrid_R, ygrid_R, wgrid, hgrid, agrid,
                          bgrid, ymedian, r, δ, λ)

# 6. Solve Working Life Problem
(v, wp) = solveWorkingLife(v, wp, wgrid, hgrid, agrid, bgrid, zgrid, stdy,
                           r, δ, λ, ρ)

# 7. Simulate Economy
(c_t, h_t, w_t, wp_t) = simulate(wp, hgrid,yit, ymedian, s_f_i, wp_R, r, λ)
