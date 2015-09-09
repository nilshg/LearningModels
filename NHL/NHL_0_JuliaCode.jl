nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import ApproXD, Grid, Distributions, Optim, PyPlot, PyCall, FastGaussQuadrature
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
include(path*"NHL/NHL_Diagnostics.jl")

# 1. Draw Income Distribution
(yit, pension) = incomeDistribution("tew207")
#(yit, pension, α, β, β_k) = incomeDistribution(agents, bs, μₐ, μᵦ, var_α,
#                                     var_β, cov_αβ, var_ɛ, var_η, ρ, y_adj, tW)

# 2. Construct individual specific belief histories
(s_f_i, stdy, k) = learning("tew207")
#(s_f_i, stdy, k) = learning(α, β, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η,
#                              var_ɛ, fpu)

# 3. Construct Grids
(xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(xpoints, apoints, bpoints, zpoints, wpoints_R, ypoints_R, power, "tew207")

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, ygrid_R, r, δ, σ)

# 5. Solve Terminal Period Problem instead of 4. and 5.
(v, wp, c_over_x) =
  solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid, zgrid, yit, r, δ)

# 6. Solve working life problem
@time (v, wp, c_over_x) =
  solveWorkingLife(v,wp, xgrid, agrid, bgrid, zgrid, stdy, k, r, δ, ρ, c_over_x)

# 7. Simulate wealth distribution
(c_t, w_t, wp_t) =
  sim(wp, xgrid, agrid, bgrid,  zgrid, wgrid_R, yit, s_f_i, pension, r, δ, σ)

agg_w = [mean(mean(w_t[:,6:15],2)) mean(mean(w_t[:,16:25],2)) mean(mean(w_t[:,26:35],2))]./mean(yit)
wpint = interpolateV(wp[:,:,:,:,39], xgrid[:,39],agrid,bgrid,zgrid)
wpex = getValue(wpint,[10.,agrid[2],bgrid[5],zgrid[3]])[1]
mrplc = mean(pension)/mean(yit)

using DataFrames
results = DataFrame(
  Metric =
    ["Parameters","xpoints","apoints", "bpoints","zpoints","xhigh40","xlow1",
     "curvature","δ",
     "Replacement rate"," ", "Results","xmax40","w_agg(26-35)/mean(y)","w_agg(36-45)/mean(y)",
     "w_agg(46-55)/mean(y)","wp[10,2.005,0.002,-0.47]"],
  Value =
    [" ",xpoints,apoints,bpoints,zpoints,round(xgrid[end,40],2),round(xgrid[1,1],2),
     power, δ, round(mrplc,2)," "," ",maximum(w_t[:,40]),round(agg_w[1],2),round(agg_w[2],2),
     round(agg_w[3],2),round(wpex,3)])

writetable("C:/Users/tew207/Desktop/results.csv", results)
println(results)
