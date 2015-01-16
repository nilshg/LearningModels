#######################################################################################
############################## PLOTS - plotting functions #############################
#######################################################################################

# Contains:
# plotyhist(agent, yit):
#     Plots the income realizations and beliefs of one agent
# plotv(v::Array{Float64, 6}, wg, hg, ag, bg, zg, ydim, a, b, z, t):
#     Plots 3D Value function with wealth on x-axis and habits, or belief about one
#     of the components of the belief vector on y-axis (ydim can be "h", "a", "b" or "z")
# plotv(v::Array{Float64, 4}, wg, hg, yg, h, y, ydim, t):
#     Plots 3D retirement Value function with habits, wealth on x-axis and habit or
#     pension income on y axis
# plotv(v::Array{Float64, 5}, wg, ag, bg, zg, ydim, t)
#     Plots 3D value function w/o habits with wealth on x-axis and belief on y-axis
# plotv(v::Array{Float64, 3}, wgrid, ygrid, t):
#     Plots 3D retirement value function w/o habits
# plotwd(w_t, periods, δ)
#     Plots

#######################################################################################

using PyCall, PyPlot, QuantEcon
@pyimport seaborn as sns

#######################################################################################

function plotv(v::Array{Float64, 6}, wg::Array, hg::Array, ag::Array, bg::Array, zg::Array,
               ydim::String, h::Int64, a::Int64, b::Int64, z::Int64, t::Int64, heading::String)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "h"
    x1g, x2g = meshgrid(wg[:, t], hg[:, t])
    x3g = v[:, :, a, b, z, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
  elseif ydim == "a"
    x1g, x2g = meshgrid(wg[:, t], ag)
    x3g = reshape(v[:, h, :, b, z, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel](L"Belief Level ($\alpha$)", fontsize=14)
  elseif ydim == "b"
    x1g, x2g = meshgrid(wg[:, t], bg)
    x3g = reshape(v[:, h, a, :, z, t]', size(v, 1), size(v, 4))
    ax[:set_ylabel](L"Belief Level ($\beta$)", fontsize=14)
  elseif ydim == "z"
    x1g, x2g = meshgrid(wg[:, t], zg)
    x3g = reshape(v[:, h, a, b, :, t]', size(v, 1), size(v, 5))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
  end

  ax[:plot_surface](x1g', x2g', x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt.show()
end

#######################################################################################

function plotv(v::Array{Float64, 4}, wg::Array, hg::Array, yg::Array, ydim::String,
               h::Int64, y::Int64, t::Int64, heading::String)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "h"
    x1g, x2g = meshgrid(wg[:, t], hg[:, t])
    x3g = v[:, :, y, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
  elseif ydim == "y"
    x1g, x2g = meshgrid(wg[:, t], yg[:, t])
    x3g = reshape(v[:, h, :, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel]("Pension Level", fontsize=14)
  end

  ax[:plot_surface](x1g, x2g, x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt.show()
end

#######################################################################################

function plotv(v::Array{Float64, 5}, wg::Array, ag::Array, bg::Array, zg::Array,
               ydim::String, a::Int64, b::Int64, z::Int64, t::Int64, heading::String)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "a"
    x1g, x2g = meshgrid(wg[:, t], ag)
    x3g = v[:, :, b, z, t]
    ax[:set_ylabel]("Belief Level (\$\\alpha\$)", fontsize=14)
  elseif ydim == "b"
    x1g, x2g = meshgrid(wg[:, t], bg)
    x3g = reshape(v[:, a, :, z, t]', size(v, 1), size(v, 3))
    ax[:set_ylabel]("Belief Level (\$\\beta\$)", fontsize=14)
  elseif ydim == "z"
    x1g, x2g = meshgrid(wg[:, t], zg)
    x3g = reshape(v[:, a, b, :, t]', size(v, 1), size(v, 4))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
  end

  ax[:plot_surface](x1g', x2g', x3g, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt.show()
end

#######################################################################################

function plotv(v::Array{Float64, 3}, wg::Array{Float64, 2}, yg::Array{Float64, 1},
               t::Int64, heading::String)

  xg, yg = meshgrid(wg[:, t], yg[:])

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  ax[:plot_surface](xg, yg, v[:, :, t]', rstride = 1, cstride = 1,
                    alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_ylabel]("Pension Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title(heading)
  plt.show()
end

#######################################################################################

function plotdistributions(w_t::Array{Float64, 2}, periods::Array, δ::Float64)

  @assert length(periods) == 4

  fig, ax = PyPlot.subplots(2, 2, sharex=true, sharey=true)
  ax[1,1][:hist](w_t[:, periods[1]], bins = 100)
  ax[1,1][:set_title]("Period "*string(periods[1]))
  ax[1,2][:hist](w_t[:, periods[2]], bins = 100)
  ax[1,2][:set_title]("Period "*string(periods[2]))
  ax[2,1][:hist](w_t[:, periods[3]], bins = 100)
  ax[2,1][:set_title]("Period "*string(periods[3]))
  ax[2,2][:hist](w_t[:, periods[4]], bins = 100)
  ax[2,2][:set_title]("Period "*string(periods[4]))
  fig[:suptitle](L"Wealth Distributions, $\delta$="*string(δ))
  plt.show()
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
  plt.show()
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
  fig[:suptitle]("Simulation History for Agent "*string(i))
  plt.legend()
  plt.show()
end

#######################################################################################

function plot2Dconfunc(t::Int64)
  fig, ax = PyPlot.subplots(2, 2)
  ax[1,1][:plot](c_over_x[:, 1, 4, 4, t], label = L"\alpha=1")
  ax[1,1][:plot](c_over_x[:, 2, 4, 4, t], label = L"\alpha=2")
  ax[1,1][:plot](c_over_x[:, 3, 4, 4, t], label = L"\alpha=3")
  ax[1,1][:plot](c_over_x[:, 4, 4, 4, t], label = L"\alpha=4")
  ax[1,1][:set_title]("c/x for "*L"\beta=4, z=4")
  ax[1,1][:legend](loc = "best")
  ax[2,1][:plot](c_over_x[:, 3, 1, 4, t], label = L"\beta=1")
  ax[2,1][:plot](c_over_x[:, 3, 4, 4, t], label = L"\beta=4")
  ax[2,1][:plot](c_over_x[:, 3, 6, 4, t], label = L"\beta=6")
  ax[2,1][:plot](c_over_x[:, 3, 9, 4, t], label = L"\beta=9")
  ax[2,1][:set_title]("c/x for "*L"\alpha=3, z=4")
  ax[2,1][:legend](loc = "best")
  ax[1,2][:plot](c_over_x[:, 3, 4, 1, t], label = L"z=1")
  ax[1,2][:plot](c_over_x[:, 3, 4, 3, t], label = L"z=2")
  ax[1,2][:plot](c_over_x[:, 3, 4, 5, t], label = L"z=3")
  ax[1,2][:plot](c_over_x[:, 3, 4, 7, t], label = L"z=4")
  ax[1,2][:set_title]("c/x for "*L"\alpha=3, \beta=4")
  ax[1,2][:legend](loc = "best")
  ax[2,2][:plot](c_over_x[:, 1, 1, 1, t], label = L"\beta=1")
  ax[2,2][:plot](c_over_x[:, 1, 4, 1, t], label = L"\beta=4")
  ax[2,2][:plot](c_over_x[:, 1, 6, 1, t], label = L"\beta=6")
  ax[2,2][:plot](c_over_x[:, 1, 9, 1, t], label = L"\beta=9")
  ax[2,2][:set_title]("c/x for "*L"\alpha=1, z=1")
  ax[2,2][:legend](loc = "best")
  fig[:suptitle]("Consumption over assets, period "*string(t))
  plt.show()
end
