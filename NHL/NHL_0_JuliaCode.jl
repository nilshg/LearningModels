using Distributed
nprocs()==Sys.CPU_THREADS || addprocs(Sys.CPU_THREADS-1)
import Distributions, FastGaussQuadrature, Interpolations, Optim, Plots, StatsBase
using SpecialFunctions
@everywhere begin
  p=("C:/Users/nils/Documents/GitHub/LearningModels/")
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
s_f_i, stdy, k = learning(α, β, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ,
                          g_t, fpu, true)

# 3. Construct Grids
xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R = grids(s_f_i, stdy, xpoints,
  apoints, bpoints, zpoints, wpoints_R, ypoints_R, power, r, g_t,
  pension, true)

# 4. Solve Retirement Problem
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ, power)

# 5. Solve Transition Period
v, wp, = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid, zgrid,
                          yit, g_t, r, δ, σ)

# 6. Solve working life problem
@time v, wp = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid, stdy, k, r, δ,
                                ρ, g_t, σ, ξ)

# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR, p)

a, b, = compare_wealth_dist(w_t,
    vec(SCF_prime_83), vec(SCF_young_83), vec(SCF_middle_83), vec(SCF_old_83))

ξ_1 = 0.005
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ, power)
# 5. Solve Transition Period
v, wp = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
# 6. Solve working life problem
@time v, wp = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, g_t, σ, ξ_1)
# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct2 = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR, p)

ξ_2 = 0.01
v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ, power)
# 5. Solve Transition Period
v, wp = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                              zgrid, yit, g_t, r, δ, σ)
# 6. Solve working life problem
@time v, wp = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, g_t, σ, ξ_2)
# 7. Simulate wealth distribution
c_t, w_t, wp_t, pct3 = sim(wp, xgrid, agrid, bgrid, zgrid, wgrid_R,
                          yit, s_f_i, pension, r, δ, σ, tR, p)

comp_statics(pct, pct2, pct3, "ξ", 0.0, ξ_1, ξ_2)
