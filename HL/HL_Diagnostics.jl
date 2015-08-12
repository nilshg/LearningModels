#######################################################################################
############################## PLOTS - plotting functions #############################
#######################################################################################

using PyPlot, PyCall
@pyimport seaborn as sns

function plotv(v::Array{Float64, 6}, wg::Array, hg::Array, ag::Array, bg::Array, zg::Array,
               ydim::String, h::Int64, a::Int64, b::Int64, z::Int64, t::Int64)

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
