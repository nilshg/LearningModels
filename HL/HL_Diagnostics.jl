#######################################################################################
############################## PLOTS - plotting functions #############################
#######################################################################################

using PyPlot, PyCall
@pyimport seaborn as sns

################################################################################

function plotv(v::Array{Float64, 6}, wg::Array, hg::Array, ag::Array, bg::Array,
  zg::Array, ydim::String, h::Int64, a::Int64, b::Int64, z::Int64, t::Int64)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "h"
    x1g, x2g = meshgrid(wg[:, t], hg[:, t])
    x3g = v[:, :, a, b, z, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
    title("v[:,:,$a,$b,$z,$t]")
  elseif ydim == "a"
    x1g, x2g = meshgrid(wg[:, t], ag)
    x3g = reshape(v[:, h, :, b, z, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel](L"Belief Level ($\alpha$)", fontsize=14)
    title("v[:,$h,:,$b,$z,$t]")
  elseif ydim == "b"
    x1g, x2g = meshgrid(wg[:, t], bg)
    x3g = reshape(v[:, h, a, :, z, t]', size(v, 1), size(v, 4))
    ax[:set_ylabel](L"Belief Level ($\beta$)", fontsize=14)
    title("v[:,$h,$a,:,$z,$t]")
  elseif ydim == "z"
    x1g, x2g = meshgrid(wg[:, t], zg)
    x3g = reshape(v[:, h, a, b, :, t]', size(v, 1), size(v, 5))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
    title("v[:,$h,$a,$b,:,$t]")
  end

  ax[:plot_surface](x1g', x2g', x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  plt[:show]()
end

################################################################################

function plotv(v::Array{Float64, 4}, wg::Array, hg::Array, yg::Array,
  ydim::String, h::Int64, y::Int64, t::Int64)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "h"
    x1g, x2g = meshgrid(wg[:, t], hg[:, t])
    x3g = v[:, :, y, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
    title("v[:,:,$y,$t]")
  elseif ydim == "y"
    x1g, x2g = meshgrid(wg[:, t], yg[:, t])
    x3g = reshape(v[:, h, :, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel]("Pension Level", fontsize=14)
    title("v[:,$h,:,$t]")
  end

  ax[:plot_surface](x1g, x2g, x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt[:show]()
end

################################################################################

function plotdistributions(w_t::Array{Float64, 2}, periods::Array, δ::Float64)

  @assert size(periods,1) == 4

  fig, axs = PyPlot.subplots(2, 2, sharex=true, sharey=true)
  for (i,ax) in enumerate(reshape(axs,4,1))
    ax[:hist](w_t[:, periods[i]], bins = 100)
    ax[:set_title]("Period "*string(periods[i]))
  end
  fig[:suptitle](L"Wealth Distributions, $\delta$="*string(δ))
  plt[:show]()
end

################################################################################

function plotdistributions(yit::Array{Float64, 2},pension::Array,periods::Array)

  @assert length(periods) == 3

  fig, ax = PyPlot.subplots(2, 2, sharex=true, sharey=true)
  ax[1,1][:hist](yit[:, periods[1]], bins = 100)
  ax[1,1][:set_title]("Period "*string(periods[1]))
  ax[1,2][:hist](yit[:, periods[2]], bins = 100)
  ax[1,2][:set_title]("Period "*string(periods[2]))
  ax[2,1][:hist](yit[:, periods[3]], bins = 100)
  ax[2,1][:set_title]("Period "*string(periods[3]))
  ax[2,2][:hist](pension, bins = 100)
  ax[2,2][:set_title]("Pension Income")
  fig[:suptitle]("Income Distributions")
  plt[:show]()
end

################################################################################

function plothistory(i::Int64, c_t::Array{Float64,2}, w_t::Array{Float64,2},
  h_t::Array{Float64,2}, yit::Array{Float64, 2}, pension::Array,
  s_f_i::Array{Float64,3}, wgrid::Array, wgrid_R::Array)

  ybelief = [exp([[1 t 1]*s_f_i[:, 1, t]][1]) for t in 1:40]

  fig, ax = PyPlot.subplots()
  ax[:plot](c_t[i, :]', label = "Consumption")
  ax[:plot](w_t[i, :]', label = "Assets")
  ax[:plot](h_t[i, :]', label = "Habits")
  ax[:plot]([yit[i, :] pension[i]*ones(30, 1)']', label = "Income")
  ax[:plot]([ybelief' pension[i]*ones(30, 1)']', label = "Belief")
  ax[:plot]([wgrid[1, :] wgrid_R[1, :]]', linestyle=":",
                                          label = "Borrowing Constraint")
  ax[:axvline](tW, linestyle = "--", color = "black")
  fig[:suptitle]("Simulation History for Agent "*string(i))
  plt.legend()
  plt[:show]()
end

################################################################################
## Simulation Results ##

function constrained_negative(w::Array{Float64,2}, wgrid::Array{Float64,2},
                              wgrid_R::Array{Float64,2})

  cnstr = zeros(Float64, size(w,2)); neg_cons = similar(cnstr)
  wmin = [wgrid[1, :] wgrid_R[1, :]]

  for t = 1:size(w,2), i = 1:size(w,1)
    abs(w[i, t] - wmin[t]) < 1e-3 ? constrained[t] += 1/1000 : 0
    c_t[i, t] < 0 ? neg_cons[t] += 1 : 0
  end
  neg_cons = neg_cons/1000.

  return constrained, neg_cons
end

################################################################################
# Variance of consumption and asset series
function crosssec_stats(c::Array{Float64,2}, w::Array{Float64,2}, y::Array,
  pension::Array, wgrid::Array, wgrid_R::Array, plot::Bool = true)

  med_c = Array(Float64, size(c,2))
  med_w = similar(med_c); mean_w = similar(med_c)
  var_c = similar(med_c); var_w = similar(var_c); var_y = similar(var_c)

  for t = 1:size(c,2)
    med_c[t] = median(c[:, t])
    med_w[t] = median(w[:, t])
    mean_w[t] = mean(w[:, t])
    var_c[t] = var(c[:, t])
    var_w[t] = var(w[:, t])
    t <= size(y,2) ? var_y[t] = var(y[:, t]) : var_y[t] = var(pension)
  end

  if plot
    bc = [wgrid[1, :] wgrid_R[1, :]]'
    fig, ax = PyPlot.subplots(2,1)
    ax[1,1][:plot](med_c, label = "Median Consumption")
    ax[1,1][:plot](med_w, label = "Median Wealth")
    ax[1,1][:plot](mean_w, label = "Mean Wealth")
    ax[1,1][:plot](bc, linestyle = ":", color = "black", label = "Borrowing Constraint")
    ax[2,1][:plot](var_c, label = "Consumption variance")
    #ax[2,1][:plot](var_w, label = "Asset variance")
    ax[2,1][:plot](var_y, label = "Income/Pension variance")
    ax[1,1][:legend](loc = "best")
    ax[2,1][:legend](loc = "best")
    ax[1,1][:set_title]("Means and Medians")
    ax[2,1][:set_title]("Variances")
    plt[:show]()
  end

  return med_c, med_w, mean_w, var_c, var_w, var_y
end

################################################################################
# Aggregate income, consumption, assets
function aggregates(c::Array{Float64,2}, w::Array{Float64,2},
  y::Array{Float64,2}, pension::Array)
  
  agg_c = Array(Float64, size(c,2))
  agg_w = similar(agg_c);  agg_y = similar(agg_c)

  for t = 1:size(c,2)
    agg_c[t] = sum(c[:, t])
    agg_w[t] = sum(w[:, t])
    t <= tW ? agg_y[t] = sum(y[:, t]) : agg_y[t] = sum(pension)
  end
  return agg_c, agg_w, agg_y
end

################################################################################
# Beliefs vs grids vs realizations
function plot_beliefs_realizations()
  ybelief = Array(Float64, (agents*bs, tW))
  ymaxbelief = Array(Float64, (tW, 1))
  yminbelief = Array(Float64, (tW, 1))
  ymingridbelief = Array(Float64, tW)
  ymaxgridbelief = Array(Float64, tW)
  yminactual = Array(Float64, tW)
  ymaxactual = Array(Float64, tW)
  for t = 1:tW
    ymingridbelief[t] = exp(agrid[1] + t*bgrid[1] + zgrid[1])
    ymaxgridbelief[t] = exp(agrid[end] + t*bgrid[end] + zgrid[end])
    yminactual[t] = minimum(yit[:, t])
    ymaxactual[t] = maximum(yit[:, t])
    for i = 1:agents*bs
      ybelief[i, t] = exp([1, t, 1]'*s_f_i[:, i, t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])
    yminbelief[t] = minimum(ybelief[:, t])
  end
  fig, ax = PyPlot.subplots(1, 2)
  ax[1,1][:plot](ymingridbelief, label = "Lowest belief on grid")
  ax[1,1][:plot](yminbelief, label = "Lowest actual belief")
  ax[1,1][:plot](yminactual, label = "Lowest actual income")
  ax[2,1][:plot](ymaxgridbelief, label = "Highest belief on grid")
  ax[2,1][:plot](ymaxbelief, label = "Highest actual belief")
  ax[2,1][:plot](ymaxactual, label = "Highest actual income")
  ax[1,1][:legend](loc="best")
  ax[2,1][:legend](loc="best")
  plt[:show]()
end
