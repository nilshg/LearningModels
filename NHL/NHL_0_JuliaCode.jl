nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import Distributions, FastGaussQuadrature, Interpolations, Optim, PyPlot, PyCall, StatsBase
@everywhere begin
  p=("C:/Users/tew207/Documents/GitHub/LearningModels/")
  include(p*"NHL/NHL_Optimizations.jl"); include(p*"NHL/NHL_Interpolations.jl")
  include(p*"Parameters.jl"); include(p*"1_Income.jl")
  include(p*"2_Learning.jl"); include(p*"NHL/NHL_3_Grid.jl")
  include(p*"NHL/NHL_4_Retirement.jl"); include(p*"NHL/NHL_5_Transition.jl")
  include(p*"NHL/NHL_6_Bellman.jl"); include(p*"NHL/NHL_7_Simulate.jl")
end
include(p*"NHL/NHL_Diagnostics.jl")

# 1. Draw Income Distribution
yit, pension, α, β, β_k, g_t, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ =
  incomeDistribution(agents, bs, tW, profile="baseline")

# 2. Construct individual specific belief histories
s_f_i, stdy, k = learning(α, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ,
                          g_t, fpu, true)

# 3. Construct Grids
xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R = grids(s_f_i, stdy, xpoints,
  apoints, bpoints, zpoints, wpoints_R, ypoints_R, power, r, tR, g_t,
  pension, true)

# 4. Solve Retirement Problem
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR)

# 5. Solve Transition Period
v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)

# 6. Solve working life problem
@time v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, ξ)

# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR)

include("C:/Users/tew207/Documents/GitHub/SCFtools/SCF_percentiles.jl")
winfriedcompare(w_t, SCF_prime_83, SCF_young_83, SCF_middle_83, SCF_old_83)

υ_2 = 0.0
v_R_2, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ_2)
v, wp, c_over_x = solveTransition(v_R_2, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
@time v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, ξ)
c_t2, w_t2, wp_t2, pct2 = sim(wp, wp_R, xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR)

υ_3 = υ
v_R_2, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ_3)
v, wp, c_over_x = solveTransition(v_R_2, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
@time v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, ξ)
c_t3, w_t3, wp_t3, pct3 = sim(wp, wp_R, xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR)

ξ_1 = 0.005
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR)
# 5. Solve Transition Period
v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
# 6. Solve working life problem
@time v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, ξ_1)
# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct2 = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR)

ξ_2 = 0.01
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR)
# 5. Solve Transition Period
v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
# 6. Solve working life problem
@time v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, ξ_2)
# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct3 = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR)

comp_statics(pct, pct2, pct3, "ξ", 0.0, ξ_1, ξ_2)
