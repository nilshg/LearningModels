include("../1_Income.jl")
include("../2_Learning.jl")
include("NHL_3_Grid.jl")
include("NHL_4_Retirement.jl")
include("NHL_5_Transition.jl")
include("NHL_7_Simulate.jl")
include("NHL_Diagnostics.jl")

nprocs() == CPU_CORES || addprocs(CPU_CORES-1)
@everywhere begin
  path=("C:/Users/tew207/My Documents/GitHub/LearningModels/")
  include(path*"Optimizations.jl")
  include(path*"Interpolations.jl")
  include(path*"Parameters.jl")
  include(path*"NHL_6_Bellman.jl")
  include(path*"../Chk_Monot.jl")
end

# 1. Draw Income Distribution
(yit, ymedian, pension) = incomeDistribution(
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning(yit, ρ, var_η, var_ɛ,
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/SNext_in.dat")

# 3. Construct Grids
(wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, apoints, bpoints, zpoints, ypoints_R, r, true,
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN code/wealth.dat",
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/wealthR.dat")

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, ygrid_R, r, δ, σ)

# 5.1 Solve Terminal Period Problem instead of 4. and 5.
(wp_1, v_1, c_over_x_1) = solveTransition(wgrid, agrid, bgrid, zgrid, r, δ, 0.0334)

@time (v_1, wp_1, c_over_x_1) =
  solveWorkingLife(v_1, wp_1, c_over_x_1, wgrid, agrid, bgrid, zgrid, stdy, k, r, δ)

(c_t_1, w_t_1, wp_t_1) =
  simulate(wp_1, wgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)
