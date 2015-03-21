nprocs()==CPU_CORES || addprocs(CPU_CORES-1)
@printf "We're using %d processors (%d workers)\n" nprocs() nworkers()

include("C:/Users/tew207/My Documents/GitHub/LearningModels/Plots.jl")
@everywhere include("C:/Users/tew207/My Documents/GitHub/LearningModels/Chk_Monot.jl")
@everywhere include("C:/Users/tew207/My Documents/GitHub/LearningModels/Optimizations.jl")
@everywhere include("C:/Users/tew207/My Documents/GitHub/LearningModels/Interpolations.jl")
@everywhere include("C:/Users/tew207/My Documents/GitHub/LearningModels/Parameters.jl")
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
φ
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
(s_f_i, stdy, k) = learning(α, β, yit, ρ, var_η, var_ɛ, guvenen_distribution)

# 3. Construct Grids
(wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R) =
  grids(s_f_i, stdy, wpoints, apoints, bpoints, zpoints,
        wpoints_R, ypoints_R, wmaxR, power, r, tR, guvenen_distribution, true)

# 4. Solve Retirement Problem
(v_R, wp_R) = solveRetirement(wgrid_R, ygrid_R, r, δ, σ)

# 5. Solve Transition Problem
#(v, wp, c_over_x) =
#  solveTransition(v_R, wgrid_R, ygrid_R, wgrid, agrid, bgrid, zgrid,
#                  yit, r, δ)

# 5.1 Solve Terminal Period Problem instead of 4. and 5.
(wp, v, c_over_x) = solveTransition(wgrid, agrid, bgrid, zgrid, r, δ)

# 6. Solve Working Life Problem
@time (v, wp, c_over_x_unc) =
  solveWorkingLife(v, wp, c_over_x, wgrid, agrid, bgrid, zgrid, stdy, k, r, δ)

# 7. Simulate consumption and wealth distributions
(c_t, w_t, wp_t) =
  simulate(wp, wgrid, agrid, bgrid,  zgrid, yit, ymedian, s_f_i, wp_R, wgrid_R,
           ygrid_R, pension, r)
