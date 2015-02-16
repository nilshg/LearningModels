include("../Plots.jl")
include("../Chk_Monot.jl")
include("../Optimizations.jl")
include("../Interpolations.jl")
include("../Parameters.jl")
include("../1_Income.jl")
include("../2_Learning.jl")
include("NHL_3_Grid.jl")
include("NHL_4_Retirement.jl")
include("NHL_5_Transition.jl")
include("NHL_6_Bellman.jl")
include("NHL_7_Simulate.jl")
include("NHL_Diagnostics.jl")

# 1. Draw Income Distribution
guvenen_distribution = true

if guvenen_distribution
  (yit, α, β, ymedian, pension) =
    incomeDistribution(
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/LaborReal.dat",
      "C:/Users/tew207/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/alfabeta.dat")
else
  (yit, α, β, ymedian, pension) =
    incomeDistribution(agents, bs, μₐ, μᵦ, var_a, var_b, var_ɛ, var_η, ρ,
                       br, tW)
end

# 2. Construct individual specific belief histories
(s_f_i, stdy) = learning(α, β, yit, ρ, var_η, var_ɛ, guvenen_distribution)

# 3. Construct Grids
(wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, stdy, wpoints, apoints, bpoints, zpoints,
        wpoints_R, ypoints_R, wmaxR, power, r, tR, guvenen_distribution, true)

# 4. Solve Retirement Problem
(v_R, wp_R, c_over_x) = solveRetirement(wgrid_R, ygrid_R, r, δ)

# 5. Solve Transition Problem
(v, wp, c_over_x) =
  solveTransition(v_R, wgrid_R, ygrid_R, wgrid, agrid, bgrid, ymedian, r, δ)

# 5.1 Solve Terminal Period Problem instead of 4. and 5.
#(wp, c, v, c_over_x) = solveTransition(wgrid, agrid, bgrid, zgrid, r, δ)

# 6. Solve Working Life Problem
(v, wp, c_over_x) =
  solveWorkingLife(v, wp, c_over_x, wgrid, agrid, bgrid, zgrid, stdy, r, δ)

# 7. Simulate consumption and wealth distributions
(c_t, w_t, wp_t) =
  simulate(wp, wgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)

(med_c, med_w, mean_w, var_c, var_w, var_y) =
  crosssec_stats(c_t, w_t, yit, pension, true)


s_f_guv_org =
  readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Guvenen FORTRAN code\\SNext_in.dat")

s_f_guv = Array(Float64, (3, 100000, 40))
for t = 1:40
  for i = 1:100000
  s_f_guv[:, i, t] = s_f_guv_org[i, 3*t-2:3*t]
  end
end

ybelief_guv = Array(Float64, (100000, 40))
error_guv = similar(ybelief_guv)
ybelief = similar(ybelief_guv)
error_my = similar(ybelief_guv)

for t = 1:40
  for i = 1:100000
    ybelief_guv[i,t] = exp([1 t 1]*s_f_guv[:, i, t])[1]
    error_guv[i,t] = yit[i,t] - ybelief_guv[i,t]
    ybelief[i,t] = exp([1 t 1]*s_f_i[:, i, t])[1]
    error_my[i,t] = yit[i,t] - ybelief[i,t]
    if t > 1
      change_guv = s_f_guv[:, i, t] - s_f_guv[:, i, t-1]
      change = s_f_i[:, i, t] - s_f_i[:, i, t-1]
    end
  end
end

fig, ax = PyPlot.subplots(2,1)
ax[1,1][:plot]([mean(abs(error_guv[:, t])) for t = 1:40], label = "Guvenen")
ax[1,1][:plot]([mean(abs(error_my[:, t])) for t = 1:40], label = "My")
ax[2,1][:hist](ybelief_guv[:, 1], bins = 100, alpha = 0.5, label = "Guvenen")
ax[2,1][:hist](ybelief[:, 1], bins = 100, alpha = 0.5, label = "My")
ax[2,1][:hist](yit[:, 1], bins = 100, alpha = 0.5, label = "True")
ax[1,1][:legend](loc = "best")
ax[2,1][:legend](loc = "best")
plt.show()
maximum(yit[:, 1])
minimum(yit[:, 1])
maximum(ybelief[:, 1])
minimum(ybelief[:, 1])
maximum(ybelief_guv[:, 1])
minimum(ybelief_guv[:, 1])
minimum(s_f_guv[1, :, 1])
minimum(s_f_guv_org[:, 1])
minimum(yit[:, 1])
