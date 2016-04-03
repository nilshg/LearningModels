nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
import Distributions, FastGaussQuadrature, Interpolations, Optim, StatsBase, NLopt
@everywhere begin
  path=("/home/tew207/")
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
  include(path*"SCF_percentiles.jl")
end

start = now()
f = open("progress.txt","a")
write(f, "\n\nSession started $start\n")
close(f)
f = open("results.txt","a")
write(f, "\n\nSession started $start\n")
close(f)

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

using StatsBase, NLopt
include(path*"NHL/NHL_Jacobian.jl")

function objective(δσ::Array{Float64,1}, grad::Array, wgrid_R=wgrid_R,
  ygrid_R=ygrid_R,r=r,xgrid=xgrid,agrid=agrid,bgrid=bgrid,zgrid=zgrid,stdy=stdy,
  k=k,s_f_i=s_f_i,pension=pension,yit=yit,data=SCF_agedetail_pctiles)
  δ = δσ[1]
  σ = δσ[2]
  υ = δσ[3]

  println("Now trying δ=$(round(δ,3)), σ=$(round(σ,3)), and υ = $(round(υ,3))")
  f = open("progress.txt","a")
  write(f, "Now trying δ=$(round(δ,3)), σ=$(round(σ,3)), and υ = $(round(υ,3))\n")
  close(f)

  v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ)

  v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid, bgrid,
                                                zgrid, yit, g_t, r, δ, σ)

  v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                      stdy, k, r, δ, ρ, c_over_x, g_t, σ, 0.0)

  c_t, w_t, wp_t, pct = sim(wp, wp_R, xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R,
                            yit, s_f_i, pension, r, δ, σ, tR)

  dif = collect(pct[10:90, 1:3] - data[10:90,2:4])
  obj = (dif'*dif)[1]
  println("Value of objective function: $obj")
  f = open("progress.txt","a")
  write(f, "Value of objective function: $obj \n")
  close(f)

  return obj
end

opt = Opt(:GN_CRS2_LM, 3)
lower_bounds!(opt, [0.91, 0.1, -2.5])
upper_bounds!(opt, [0.99, 3.0, -0.1])
min_objective!(opt, objective)
ftol_rel!(opt, 0.001)
init = [0.956 + 0.01*randn(), 2.0 + rand(), -1.5 + rand()]
(minf, minx, ret) = NLopt.optimize(opt, init)
f = open("progress.txt","a")
write(f,"Found minimum for baseline profile! Calculating standard errors...\n")
close(f)
stderr = jacobian(minx[1], minx[2], minx[3])
println("Found minimum for baseline profile at $minf, δ is $(minx[1]) ($(stderr[1])), σ is $(minx[2]) ($(stderr[2])), υ is $(minx[3]) ($(stderr[3]))")
f = open("progress.txt","a")
write(f, "Found minimum without profile at $minf, δ is $(minx[1]) ($(stderr[1])), σ is $(minx[2]) ($(stderr[2])), υ is $(minx[3]) ($(stderr[3]))\n")
close(f)

for prof in ["psid","bhps_grosslab","bhps_netlab", "bhps_netinc"]
  f = open("progress.txt","a")
  write(f, "Now calculating for $prof\n")
  close(f)
  f = open("results.txt","a")
  write(f, "Now calculating for $prof\n")
  close(f)
  yit, pension, α, β, β_k, g_t, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ =
    incomeDistribution(agents, bs, tW, profile=prof)
  s_f_i, stdy, k = learning(α, β_k, yit, ρ, var_α, var_β, cov_αβ, var_η, var_ɛ, g_t, fpu)
  xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R = grids(s_f_i, stdy, xpoints,
    apoints, bpoints, zpoints, wpoints_R, ypoints_R, power, r, tR, g_t,
    pension, true)
  opt = Opt(:GN_CRS2_LM, 3)
  lower_bounds!(opt, [0.91, 0.1, -2.5])
  upper_bounds!(opt, [0.99, 3.0, -0.1])
  min_objective!(opt, objective)
  ftol_rel!(opt, 0.001)
  init = [0.956 + 0.01*randn(), 2.0 + rand(), -1.5 + rand()]
  (minf, minx, ret) = NLopt.optimize(opt, init)
  f = open("progress.txt","a")
  write(f,"Found minimum for $prof profile! Calculating standard errors...\n")
  close(f)
  stderr = jacobian(minx[1], minx[2], minx[3])
  println("Found minimum for $prof profile at $minf, δ is $(round(minx[1],5)) ($(round(stderr[1],5))), σ is $(round(minx[2],3)) ($(round(stderr[2],5))), υ is $(round(minx[3],3)) ($(round(stderr[3],3)))")
  f = open("results.txt","a")
  write(f, "Found minimum for $prof profile at $minf, δ is $(minx[1]) ($(stderr[1])), σ is $(minx[2]) ($(stderr[2])), υ is $(minx[3]) ($(stderr[3]))\n")
  close(f)
end

start = now()
f = open("progress.txt","a")
write(f, "\nSession ended $start\n")
write(f, repeat("*",80))
close(f)
f = open("results.txt","a")
write(f, "\nSession started $start\n")
write(f, repeat("*",80))
close(f)
