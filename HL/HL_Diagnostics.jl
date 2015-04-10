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
