path=("C:/Users/tew207/My Documents/GitHub/LearningModels/")
include(path*"Optimizations.jl")
include(path*"Interpolations.jl")
include(path*"Parameters.jl")
include(path*"NHL/NHL_6_Bellman.jl")
include(path*"Chk_Monot.jl")
include(path*"1_Income.jl")
include(path*"2_Learning.jl")
include(path*"NHL/NHL_3_Grid.jl")
include(path*"NHL/NHL_4_Retirement.jl")
include(path*"NHL/NHL_5_Transition.jl")
include(path*"NHL/NHL_7_Simulate.jl")
include(path*"NHL/NHL_Diagnostics.jl")

# 1. Draw Income Distribution
(yit, ymedian, pension) = incomeDistribution(
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
  "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning(
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code")

# 3. Construct Grids
(wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, apoints, bpoints, zpoints, ypoints_R, r, true,
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN code/wealth.dat",
    "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/wealthR.dat")

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, convert(Array,ygrid_R), r, δ, σ)

# 5. Solve Terminal Period Problem instead of 4. and 5.
(wp, v, c_over_x) =
  solveTransition(v_R, wgrid_R, ygrid_R, wgrid, agrid, bgrid,
  zgrid, yit, r, δ)

# 6. Solve working life problem
(v, wp, c_over_x) =
  solveWorkingLife(v, wp, c_over_x, wgrid, agrid, bgrid, zgrid, stdy, k, r, δ)

# 7. Simulate wealth distribution
(c_t, w_t, wp_t) =
  simulate(wp, wgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)
