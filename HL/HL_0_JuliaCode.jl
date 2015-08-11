using DataFrames, PyPlot
include("../Plots.jl")
include("../1_Income.jl")
include("../2_Learning.jl")
include("HL_3_Grid.jl")
include("HL_4_Retirement.jl")
include("HL_5_Transition.jl")
include("HL_7_Simulate.jl")

nprocs() == CPU_CORES || addprocs(CPU_CORES-1)
@everywhere begin
  using Distributions, Grid, Optim, QuantEcon
  path=("C:/Users/tew207/Documents/GitHub/LearningModels/")
  include(path*"Optimizations.jl")
  include(path*"Interpolations.jl")
  include(path*"Parameters.jl")
  include(path*"Chk_Monot.jl")
  include(path*"/HL/HL_6_Bellman.jl")
end

# 1. Draw Income Distribution
(yit, ymedian, pension) = incomeDistribution(
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning(yit, ρ, var_η, var_ɛ,
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/SNext_in.dat")

# 3. Construct Grids
(wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R) =
    grids(s_f_i, stdy, wpoints, hpoints, apoints, bpoints, zpoints,
           wpoints_R, hpoints_R, ypoints_R, wmaxR, power, r, tR, true, true)

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, hgrid_R, ygrid_R, r, δ, λ)

# 5. Solve Transition Problem
(v, wp) = solveTransition(v_R, wgrid_R, hgrid_R, ygrid_R, wgrid, hgrid, agrid, bgrid,
                          zgrid, ymedian, r, δ, λ)

# 6. Solve Working Life Problem
(v, wp) = solveWorkingLife(v, wp, wgrid, hgrid, agrid, bgrid, zgrid,
                           stdy, r, k, δ, λ, ρ)

# 7. Simulate Economy
(c_t, h_t, w_t, wp_t) = simulate(wp, wgrid, hgrid, agrid, bgrid, zgrid,
                                 wgrid_R, hgrid_R, ygrid_R, yit, ymedian,
                                 s_f_i, wp_R, r, λ)
