nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import ApproXD, Distributions, Optim, PyPlot, PyCall, FastGaussQuadrature
@everywhere begin
  path=("C:/Users/Nils-Holger/Documents/GitHub/LearningModels/")
  include(path*"NHL/NHL_Optimizations.jl")
  include(path*"NHL/NHL_Interpolations.jl")
  include(path*"Parameters.jl")
  include(path*"1_Income.jl")
  include(path*"2_Learning.jl")
  include(path*"NHL/NHL_3_Grid.jl")
  include(path*"NHL/NHL_4_Retirement.jl")
  include(path*"NHL/NHL_5_Transition.jl")
  include(path*"NHL/NHL_6_Bellman.jl")
  include(path*"NHL/NHL_7_Simulate.jl")
end
include(path*"NHL/NHL_Diagnostics.jl")

# Life cycle profile from PSID 1968-1996
g_t = [8.36 + 0.07*t - 0.15*t^2/100 for t = 1:40]

# 1. Draw Income Distribution
#(yit, pension) = incomeDistribution("tew207")
yit, pension, α, β, β_k = incomeDistribution(agents, bs, μₐ, μᵦ, var_α,
                                     var_β, cov_αβ, var_ɛ, var_η, ρ, g_t, tW)

# 2. Construct individual specific belief histories
#(s_f_i, stdy, k) = learning("tew207")
s_f_i, stdy, k = learning(α, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ,
                              g_t, fpu)

# 3. Construct Grids
xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R = grids(s_f_i, stdy, xpoints,
  apoints, bpoints, zpoints, wpoints_R, ypoints_R, power, r, tR, g_t,
  pension, true)

# 4. Solve Retirement Problem
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR)

# 5. Solve Terminal Period Problem instead of 4. and 5.
v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid,
                                    bgrid, zgrid, yit, g_t, r, δ)

# 6. Solve working life problem
v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                            stdy, k, r, δ, ρ, c_over_x, g_t)

# 7. Simulate wealth distribution
c_t, w_t, wp_t = sim(wp, xgrid, agrid, bgrid,  zgrid, wgrid_R, yit, s_f_i,
                                                        pension, r, δ, σ, tR)
