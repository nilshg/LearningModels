#######################################################################################
############################## PLOTS - plotting functions #############################
#######################################################################################

# Contains:
# plotYhist(agent, Yit):
#     Plots the income realizations and beliefs of one agent
# plotV(V::Array{Float64, 6}, wg, hg, ag, bg, zg, ydim, a, b, z, t):
#     Plots 3D value function with wealth on x-axis and habits, or belief about one
#     of the components of the belief vector on y-axis (ydim can be "h", "a", "b" or "z")
# plotV(V::Array{Float64, 4}, wg, hg, yg, h, y, ydim, t):
#     Plots 3D retirement value function with habits, wealth on x-axis and habit or 
#     pension income on y axis
# plotV(V::Array{Float64, 5}, wg, ag, bg, zg, ydim, t)
#     Plots 3D value function w/o habits with wealth on x-axis and belief on y-axis
# plotV(V::Array{Float64, 3}, wgrid, ygrid, t):
#     Plots 3D retirement value function w/o habits

#######################################################################################

function plotYhist(agent::Int64, Yit::Array{Float64, 2})
  ybelief = zeros(1,40)
  yupper = zeros(1,40)
  ylower = zeros(1,40)

  for t in [1:40]
      a = H[:, t]' * S_f_i[:, agent, t]
      ybelief[1, t] = a[1]
      yupper[1, t] = a[1] + 3*stdy[t]
      ylower[1, t] = a[1] - 3*stdy[t]
  end
  fig, aPyPlot.subplots()
  ax[:plot]([1:40], exp(ybelief'), label = L"$\hat{y}^i_{t|t-1}$")
  ax[:plot]([1:40], exp(yupper'), ":", color = "cornflowerblue")
  ax[:plot]([1:40], exp(ylower'), ":", color = "cornflowerblue")
  ax[:plot]([1:40], Yit[agent, :]', label = L"$y^i_t$")
  ax[:plot]([1:40], exp(ybelief') - Yit[agent, :]', label = L"$\xi^i_t$")
  title("Income Realizations vs. Beliefs, agent " * string(agent))
  ax[:legend](loc=0)
end

#######################################################################################

function plotV(V::Array{Float64, 6}, wg::Array, hg::Array, ag::Array, bg::Array, zg:Array, 
               ydim::String, h::Int64, a::Int64, b::Int64, z::Int64, t::Int64)

  if ydim == "h"
    xg, yg = meshgrid(wg[:, t], hg[:, t])
    zg = V[:, :, a, b, z, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
  elseif ydim == "a"
    xg, yg = meshgrid(wg[:, t], ag[:, t])
    zg = reshape(V[:, h, :, b, z, t]', size(V, 1), size(V, 3))
    ax[:set_ylabel](L"Belief Level ($α$)", fontsize=14)
  elseif ydim == "b"
    xg, yg = meshgrid(wg[:, t], bg[:, t])
    zg = reshape(V[:, h, a, :, z, t]', size(V, 1), size(V, 4))
    ax[:set_ylabel](L"Belief Level ($β$)", fontsize=14)
  elseif ydim == "z"
    xg, yg = meshgrid(wg[:, t], zg[:, t])
    zg = reshape(V[:, h, a, b, :, t]', size(V, 1), size(V, 5))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
  end
  
  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  ax[:plot_surface](xg, yg, zg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title("Value Function for period " * string(t))
  plt.show()
end

#######################################################################################

function plotV(V::Array{Float64, 4}, wg::Array, hg::Array, yg::Array,
               ydim::String, h, y, t::Real)

  if ydim == "h"
    xg, yg = meshgrid(wg[:, t], hg[:, t])
    zg = V[:, :, y, t]'
    ax[:set_ylabel]("Habit Level", fontsize=14)
  elseif ydim == "y"
    xg, yg = meshgrid(wg[:, t], hg[:, t])
    zg = reshape(V[:, h, :, t]', size(V, 1), size(V, 3))
    ax[:set_ylabel]("Pension Level", fontsize=14)
  end    

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  ax[:plot_surface](xg, yg, zg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title("Value Function for retirement period " * string(t))
  plt.show()
end

#######################################################################################

function plotV(V::Array{Float64, 5}, wg::Array, ag::Array, bg::Array, zg::Array,
               ydim::String, a::Int64, b::Int64, z::Int64, t::Int64)

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  if ydim == "a"
    xg, yg = meshgrid(wgrid[:, t], agrid[:, t])
    zg = V[:, :, b, z, t]
    ax[:set_ylabel](L"Belief Level (\$\\alpha\$)", fontsize=14)
  elseif ydim == "b"
    xg, yg = meshgrid(wgrid[:, t], bgrid[:, t])
    zg = reshape(V[:, a, :, z, t]', size(V, 1), size(V, 3))
    ax[:set_ylabel](L"Belief Level (\$\\beta\$)", fontsize=14)
  elseif ydim == "z"
    xg, yg = meshgrid(wgrid[:, t], zgrid[:, t])
    zg = reshape(V[:, a, b, :, t]', size(V, 1), size(V, 4))
    ax[:set_ylabel]("Belief Level (z)", fontsize=14)
  end

  ax[:plot_surface](xg', yg', zg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title("Value Function for period " * string(t))
  plt.show()
end

#######################################################################################

function plotV(V::Array{Float64, 3}, wg::Array, yg::Array, t::Int64)

  xg, yg = meshgrid(wg[:, t], yg[:])

  fig = figure(figsize=(10,8))
  ax = fig[:add_subplot](111, projection="3d")

  ax[:plot_surface](xg, yg, V[:, :, t]', rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
  ax[:set_xlabel]("Wealth Level", fontsize=14)
  ax[:set_ylabel]("Pension Level", fontsize=14)
  ax[:set_zlabel]("Value", fontsize=14)
  title("Value Function for retirement period " * string(t))
  plt.show()
end

#######################################################################################

function plotV2d(V::Array{Float64, 6}, wgrid::Array, hgrid::Array, agrid::Array,
                 bgrid::Array, zgrid::Array, w::Real, h::Real, y::Real, a::Real, b::Real,
                 z::Real, t::Real; dim="w")

  fig = figure(figsize=(10,8))
  fig[:add_subplot](111)
  if dim == "w"
    plot(wgrid[:, t], V[:, h, a, b, z, t]')
    ax[:set_xlabel]("Wealth Level")
    ax[:set_ylabel]("Value")
    title("Value Function for period " * string(t))
  elseif dim == "h"
    plot(hgrid[:, t], V[w, :, a, b, z, t]')
    ax[:set_xlabel]("Habit Level")
    ax[:set_ylabel]("Value")
    title("Value Function for period " * string(t))
  elseif dim == "a"
    plot(agrid[:, t], V[w, h, :, b, z, t]')
    ax[:set_xlabel]("Belief Level (α)")
    ax[:set_ylabel]("Value")
    title("Value Function for period " * string(t))
  elseif dim == "b"
    plot(bgrid[:, t], V[w, h, a, :, z, t]')
    ax[:set_xlabel]("Belief Level (β)")
    ax[:set_ylabel]("Value")
    title("Value Function for period " * string(t))
  elseif dim == "b"
    plot(zgrid[:, t], V[w, h, a, b, :, t]')
    ax[:set_xlabel]("Belief Level (z)")
    ax[:set_ylabel]("Value")
    title("Value Function for period " * string(t))
  end
  plt.show()
end

#######################################################################################

function plotV2d(V::Array{Float64, 4}, wgrid::Array, hgrid::Array, ygrid::Array,
                 w::Real, h::Real, y::Real, t::Real; dim="w")

  fig = figure(figsize=(10,8))
  fig[:add_subplot](111)
  if dim == "w"
    plot(wgrid[:, t], V[:, h, y, t]')
    ax[:set_xlabel]("Wealth Level")
    ax[:set_ylabel]("Value")
    title("Retirement Value Function for period " * string(t))
  elseif dim == "h"
    plot(hgrid[:, t], V[w, :, a, b, z, t]')
    ax[:set_xlabel]("Habit Level")
    ax[:set_ylabel]("Value")
    title("Retirement Value Function for period " * string(t))
  elseif dim == "a"
    plot(agrid[:, t], V[w, h, :, b, z, t]')
    ax[:set_xlabel]("Belief Level (α)")
    ax[:set_ylabel]("Value")
    title("Retirement Value Function for period " * string(t))
  elseif dim == "b"
    plot(bgrid[:, t], V[w, h, a, :, z, t]')
    ax[:set_xlabel]("Belief Level (β)")
    ax[:set_ylabel]("Value")
    title("Retirement Value Function for period " * string(t))
  elseif dim == "b"
    plot(zgrid[:, t], V[w, h, a, b, :, t]')
    ax[:set_xlabel]("Belief Level (z)")
    ax[:set_ylabel]("Value")
    title("Retirement Value Function for period " * string(t))
  end
  plt.show()
end
