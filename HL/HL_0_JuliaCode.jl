nprocs() == CPU_CORES || addprocs(CPU_CORES-1)
import ApproXD, Grid, Distributions, Optim, QuantEcon, PyPlot, PyCall
@everywhere begin
  path=("C:/Users/tew207/Documents/GitHub/LearningModels/")
  include(path*"Optimizations.jl")
  include(path*"Interpolations.jl")
  include(path*"Parameters.jl")
  include(path*"1_Income.jl")
  include(path*"2_Learning.jl")
  include(path*"HL/HL_3_Grid.jl")
  include(path*"HL/HL_4_Retirement.jl")
  include(path*"HL/HL_5_Transition.jl")
  include(path*"HL/HL_6_Bellman.jl")
  include(path*"HL/HL_7_Simulate.jl")
  include(path*"Chk_Monot.jl")
end

# 1. Draw Income Distribution
(yit, ymedian, pension) = incomeDistribution("tew207")

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning("tew207")

# 3. Construct Grids
(xgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R) =
    grids(s_f_i, stdy, wpoints, hpoints, apoints, bpoints, zpoints,
           wpoints_R, hpoints_R, ypoints_R, wmaxR, power, r, tR, true)

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, hgrid_R, ygrid_R, r, δ, λ)

# 5. Solve Transition Problem
(v, wp, c_over_x) =
  solveTransition(v_R, wgrid_R, hgrid_R, ygrid_R, xgrid, hgrid, agrid, bgrid,
                  zgrid, ymedian, r, δ, λ)

# 6. Solve Working Life Problem
(v, wp, c_over_x) = solveWorkingLife(v, wp, xgrid, hgrid, agrid, bgrid, zgrid,
                                     stdy, k, r, δ, λ, ρ, c_over_x)

# 7. Simulate Economy
(c_t, h_t, w_t, wp_t) = sim(wp, xgrid, hgrid, agrid, bgrid, zgrid,
                  wgrid_R, hgrid_R, ygrid_R, yit, s_f_i, wp_R, pension, r, λ)
