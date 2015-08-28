################################################################################
################################  DIAGNOSTICS  #################################
################################################################################

using PyCall, PyPlot, QuantEcon
@pyimport seaborn as sns

#######################################################################################

function plotv(v::Array{Float64, 5}, wg::Array, ag::Array, bg::Array, zg::Array,
               ydim::String, a::Int64, b::Int64, z::Int64, t::Int64)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "a"
    x1g, x2g = meshgrid(wg[:, t], ag)
    x3g = v[:, :, b, z, t]
    ax[:set_ylabel]("Belief Level (\$\\alpha\$)", fontsize=14)
    title("v[:,:,$b,$z,$t]")
  elseif ydim == "b"
    x1g, x2g = meshgrid(wg[:, t], bg)
    x3g = reshape(v[:, a, :, z, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel]("Belief Level (\$\\beta\$)", fontsize=14)
    title("v[:,$a,:,$z,$t]")
  elseif ydim == "z"
    x1g, x2g = meshgrid(wg[:, t], zg)
    x3g = reshape(v[:, a, b, :, t]', size(v, 1), size(v, 4))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
    title("v[:,$a,$b,:,$t]")
  end

  ax[:plot_surface](x1g', x2g', x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  plt[:show]()
end

#######################################################################################

function plotv(v::Array{Float64, 3}, wg::Array{Float64, 2}, yg::Array{Float64, 1},
               t::Int64; heading::String = "")

  xg, yg = meshgrid(wg[:, t], yg[:])

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  ax[:plot_surface](xg, yg, v[:, :, t]', rstride = 1, cstride = 1,
                    alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_ylabel]("Pension Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt[:show]()
end

#######################################################################################

# Plot value and policy functions
function plot_value_policy(tw::Int64, tr::Int64)
  plotv(v, wgrid, agrid, bgrid, zgrid, "b", 3, 7, 1, tw,
      "Value Function, period "*string(tw))
  plotv(wp, wgrid, agrid, bgrid, zgrid, "b", 3, 5, 4, tw,
      "Policy Function, period "*string(tw))
  plotv(v_R, wgrid_R, ygrid_R, tr,
      "Retirement value function, period "*string(tr))
  plotv(wp_R, wgrid_R, ygrid_R, tr,
      "Retirement policy function, period "*string(tr))
end

#######################################################################################

function plotdistributions(w_t::Array{Float64, 2}, periods::Array, δ::Float64)

  @assert size(periods,1) == 4

  fig, axes = PyPlot.subplots(2, 2, sharex=true, sharey=true)
  for (i,ax) in enumerate(reshape(axes,4,1))
    ax[:hist](w_t[:, periods[i]], bins = 100)
    ax[:set_title]("Period "*string(periods[i]))
  end
  fig[:suptitle](L"Wealth Distributions, $\delta$="*string(δ))
  plt[:show]()
end

#######################################################################################

function plotdistributions(yit::Array{Float64, 2}, pension::Array, periods::Array)

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

#######################################################################################

function plothistory(i::Int64, c_t::Array{Float64, 2}, w_t::Array{Float64, 2},
                     yit::Array{Float64, 2}, pension::Array, s_f_i::Array{Float64, 3},
                     wgrid::Array, wgrid_R::Array, tW::Int64, tR::Int64)

  ybelief = [exp([[1 t 1]*s_f_i[:, 1, t]][1]) for t in 1:40]
  fig, ax = PyPlot.subplots()
  ax[:plot](c_t[i, :]', label = "Consumption")
  ax[:plot](w_t[i, :]', label = "Assets")
  ax[:plot]([yit[i, :] pension[i]*ones(tR, 1)']', label = "Income")
  ax[:plot]([ybelief' pension[i]*ones(tR, 1)']', label = "Belief")
  ax[:plot]([wgrid[1, :] wgrid_R[1, :]]', linestyle=":", label = "Borrowing Constraint")
  ax[:axvline](tW, linestyle = "--", color = "black")
  ax[:legend](loc="best")
  fig[:suptitle]("Simulation History for Agent "*string(i))
  plt[:show]()
end

#######################################################################################
## Simulation Results ##

function constrained_negative(w::Array{Float64,2}, wgrid::Array{Float64,2},
                              wgrid_R::Array{Float64,2})

  constrained = zeros(Float64, size(w,2))
  wmin = [wgrid[1, :] wgrid_R[1, :]]

  for t = 1:size(wgrid,2), i = 1:size(wp,1)
    if abs(w[i, t] - wmin[t]) < 1e-3
      constrained[t] += 1/1000
    end
  end

  neg_cons = similar(constrained)
  for t = 1:70, i = 1:size(w,1)
    if c_t[i, t] < 0
      neg_cons[t] += 1
    end
  end
  neg_cons = neg_cons/1000

  return constrained, neg_cons
end
################################################################################
# Variance of consumption and asset series
function crosssec_stats(c::Array{Float64,2}, w::Array{Float64,2}, y::Array,
                        pension::Array, wgrid::Array, wgrid_R::Array,
                        plot::Bool = true)

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
    ax[2,1][:plot](var_w, label = "Asset variance")
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

################################################################################
function plot2Dconfunc(c_x::Array{Float64, 5}, t::Int64, wgrid::Array, a::Int64,
  b::Int64; figtext::String = "")

  α_mid = convert(Int64, round(size(c_x,2)/2))
  β_mid = convert(Int64, round(size(c_x,3)/2))
  z_mid = convert(Int64, round(size(c_x,4)/2))
  α_hi = size(c_x,2)
  β_hi = size(c_x,3)
  z_hi = size(c_x,4)

  fig, ax = PyPlot.subplots(2, 2)
  ax[1,1][:plot](wgrid[a:b,t], c_x[a:b, 1, β_mid, z_mid, t], label = "α=1")
  ax[1,1][:plot](wgrid[a:b,t], c_x[a:b, α_mid, β_mid, z_mid, t], label = "α=2")
  ax[1,1][:plot](wgrid[a:b,t], c_x[a:b,end,β_mid,z_mid,t], label="α="*string(α_hi))
  #ax[1,1][:plot](wgrid[a:b,t],wgrid[a:b,t], linestyle=":", color="black")
  ax[1,1][:set_title]("β="*string(β_mid)*", z="*string(z_mid))
  ax[1,1][:legend](loc = "best")
  ax[2,1][:plot](wgrid[a:b,t], c_x[a:b, α_mid, 1, 4, t], label = "β=1")
  ax[2,1][:plot](wgrid[a:b,t], c_x[a:b, α_mid, 4, 4, t], label = "β=4")
  ax[2,1][:plot](wgrid[a:b,t], c_x[a:b, α_mid, 6, 4, t], label = "β=6")
  ax[2,1][:plot](wgrid[a:b,t], c_x[a:b,α_mid,end,4,t],label = "β="*string(β_hi))
  #ax[2,1][:plot](wgrid[a:b,t],wgrid[a:b,t], linestyle=":", color="black")
  ax[2,1][:set_title]("α="*string(α_mid)*", z="*string(z_mid))
  ax[2,1][:legend](loc = "best")
  ax[1,2][:plot](wgrid[a:b,t], c_x[a:b, α_mid, β_mid, 1, t], label = "z=1")
  ax[1,2][:plot](wgrid[a:b,t], c_x[a:b, α_mid, β_mid, 3, t], label = "z=2")
  ax[1,2][:plot](wgrid[a:b,t], c_x[a:b, α_mid, β_mid, 5, t], label = "z=3")
  ax[1,2][:plot](wgrid[a:b,t], c_x[a:b,α_mid,β_mid,end,t],label="z="*string(z_hi))
  #ax[1,2][:plot](wgrid[a:b,t],wgrid[a:b,t], linestyle=":", color="black")
  ax[1,2][:set_title]("α="*string(α_mid)*", "*"β="*string(β_mid))
  ax[1,2][:legend](loc = "best")
  ax[2,2][:plot](wgrid[a:b,t], c_x[a:b, 1, 1, 1, t], label = "β=1")
  ax[2,2][:plot](wgrid[a:b,t], c_x[a:b, 1, 4, 1, t], label = "β=4")
  ax[2,2][:plot](wgrid[a:b,t], c_x[a:b, 1, 6, 1, t], label = "β=6")
  ax[2,2][:plot](wgrid[a:b,t], c_x[a:b, 1, 8, 1, t], label = "β="*string(β_hi))
  #ax[2,2][:plot](wgrid[a:b,t],wgrid[a:b,t], linestyle=":", color="black")
  ax[2,2][:set_title]("α=1, z=1")
  ax[2,2][:legend](loc = "best")
  fig[:suptitle]("c/x, period "*string(t))
  fig[:text](0.01, 0.01, figtext)
  plt[:show]()
end
