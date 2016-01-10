nprocs() == CPU_CORES || addprocs(CPU_CORES-1)
import Distributions, FastGaussQuadrature, Interpolations, Optim, PyPlot, PyCall, StatsBase
@everywhere begin
  p=("C:/Users/tew207/Documents/GitHub/LearningModels/")
  include(p*"HL/HL_Interpolations.jl"); include(p*"HL/HL_Optimizations.jl");
  include(p*"Parameters.jl"); include(p*"1_Income.jl")
  include(p*"2_Learning.jl"); include(p*"HL/HL_3_Grid.jl")
  include(p*"HL/HL_4_Retirement.jl"); include(p*"HL/HL_5_Transition.jl")
  include(p*"HL/HL_6_Bellman.jl"); include(p*"HL/HL_7_Simulate.jl")
  include(p*"Chk_Monot.jl")
end

# 1. Draw Income Distribution
yit, pension, α, β, β_k, g_t, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ =
  incomeDistribution(agents, bs, tW, profile="baseline")

# 2. Construct individual specific belief histories
s_f_i, stdy, k = learning(α, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ,
                          g_t, fpu, true)

# 3. Construct Grids
(xgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R) =
    grids(s_f_i, stdy, xpoints, hpoints, apoints, bpoints, zpoints,
           wpoints_R, hpoints_R, ypoints_R, wmaxR, power, r, tR, true)

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, hgrid_R, ygrid_R, r, δ, λ)

# 5. Solve Transition Problem
(v, wp, c_over_x) =
  solveTransition(v_R, wgrid_R, hgrid_R, ygrid_R, xgrid, hgrid, agrid, bgrid,
                  zgrid, r, δ, λ)

# 6. Solve Working Life Problem
(v, wp, c_over_x) = solveWorkingLife(v, wp, xgrid, hgrid, agrid, bgrid, zgrid,
                                     stdy, k, r, δ, λ, ρ, c_over_x)

# 7. Simulate Economy
(c_t, h_t, w_t, wp_t) = sim(wp, xgrid, hgrid, agrid, bgrid, zgrid,
                  wgrid_R, hgrid_R, ygrid_R, yit, s_f_i, wp_R, pension, r, λ)

include("C:/Users/tew207/Documents/GitHub/SCFtools/SCF_percentiles.jl")
winfriedcompare(w_t./mean(yit), SCF_prime_83, SCF_young_83, SCF_middle_83, SCF_old_83)
