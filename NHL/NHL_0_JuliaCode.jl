nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import ApproXD, Grid, Distributions, Optim, QuantEcon, PyPlot, PyCall
@everywhere begin
  path=("C:/Users/tew207/My Documents/GitHub/LearningModels/")
  include(path*"Optimizations.jl")
  include(path*"Interpolations.jl")
  include(path*"Parameters.jl")
  include(path*"1_Income.jl")
  include(path*"2_Learning.jl")
  include(path*"NHL/NHL_3_Grid.jl")
  include(path*"NHL/NHL_4_Retirement.jl")
  include(path*"NHL/NHL_5_Transition.jl")
  include(path*"NHL/NHL_6_Bellman.jl")
  include(path*"NHL/NHL_7_Simulate.jl")
end

# 1. Draw Income Distribution
(yit, ymedian, pension) = incomeDistribution("tew207")

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning("tew207")

# 3. Construct Grids
(xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, apoints, bpoints, zpoints, ypoints_R, r, true, "tew207")

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, ygrid_R, r, δ, σ)

# 5. Solve Terminal Period Problem instead of 4. and 5.
(v, wp, c_over_x) =
  solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid, zgrid, yit, r, δ)

# 6. Solve working life problem
@time (v, wp, c_over_x) =
  solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid, stdy, k, r, δ, c_over_x)

# 7. Simulate wealth distribution
(c_t, w_t, wp_t) =
  sim(wp, xgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)
