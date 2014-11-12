######################################################################
using PyPlot
using Grid
using Dierckx
#using ApproXD

meshgrid(x::Vector, y::Vector) = (repmat(x, 1, length(y))',
                                  repmat(y, 1, length(x)))
######################################################################

f(x::Float64) = -(x^(-1.0))

zmin = -0.1
zmax = 1000.0
ymin = 0.2
ymax = 30.0
xmin = ymin + zmin 
xmax = ymax + zmax

xpoints = 100
ypoints = 20
zpoints = 15

xgrid = linspace(xmin, xmax, xpoints)
ygrid = linspace(ymin, ymax, ypoints)
zgrid = linspace(zmin, zmax, zpoints)

######################################################################

V_2D = zeros(ypoints, zpoints, 1)
V_1D = zeros(xpoints, 1)

# Values on grid, 1 dimension
for x = 1:xpoints
  xt = xgrid[x]
  V_1D[x] = f(xt)
end

# Values on grid, 2 dimensions
for z = 1:zpoints
  for y = 1:ypoints
    zt = zgrid[z]
    yt = ygrid[y]
    V_2D[y, z, 1] = f(yt+zt)
  end
end

######################################################################

# Interpolate in 1 dimension
xrange = range(xgrid[1], xgrid[2]-xgrid[1], xpoints)

# Grid
V_1D_int_Lin = CoordInterpGrid(xrange, V_1D[:, 1], BCnearest, InterpLinear)
V_1D_int_Quad = CoordInterpGrid(xrange, V_1D[:, 1], BCnearest, InterpQuadratic)
#V_1D_int_Cub = CoordInterpGrid(xrange, V_1D[:, 1], BCnil, InterpCubic)

# ApproXD
# V_1D_int_ApproXD = 

# Dierckx
V_1D_int_Dierckx = Spline1D(xgrid, V_1D[:, 1])

######################################################################  
  
# Interpolate in 2 dimensions
yrange = range(ygrid[1], ygrid[2]-ygrid[1], ypoints)
zrange = range(zgrid[1], zgrid[2]-zgrid[1], zpoints)

# Grid
V_2D_int_Lin = CoordInterpGrid((yrange, zrange), V_2D[:, :, 1], BCnearest, InterpLinear)
V_2D_int_Quad = CoordInterpGrid((yrange, zrange), V_2D[:, :, 1], BCnearest, InterpQuadratic)
#V_2D_int_Cub = CoordInterpGrid((yrange, zrange), V_2D[:, :, 1], BCnan, InterpCubic)

# ApproXD
# V_2D_int_ApproXD = 

# Dierckx
V_2D_int_Dierckx = Spline2D(ygrid, zgrid, V_2D[:, :, 1])

######################################################################

# Evaluate Interpolants Off-Grid
xoffgrid = linspace(xgrid[1], xgrid[end], 5*xpoints)
yoffgrid = linspace(ygrid[1], ygrid[end], 5*ypoints)
zoffgrid = linspace(zgrid[1], zgrid[end], 5*zpoints)

V_1D_offgrid_Lin = zeros(5*xpoints, 1)
V_1D_offgrid_Quad = zeros(5*xpoints, 1)
V_1D_offgrid_Dierckx = zeros(5*xpoints, 1)
V_1D_actual = zeros(5*xpoints, 1)

V_2D_offgrid_Lin =  zeros(5*ypoints, 5*zpoints, 1)
V_2D_offgrid_Quad = zeros(5*ypoints, 5*zpoints, 1)
V_2D_offgrid_Dierckx = zeros(5*ypoints, 5*zpoints, 1)
V_2D_actual = zeros(5*ypoints, 5*zpoints, 1)

for x = 1:5*xpoints
  xt = xoffgrid[x]
  V_1D_offgrid_Lin[x] = V_1D_int_Lin[xt]
  V_1D_offgrid_Quad[x] = V_1D_int_Quad[xt]
  V_1D_offgrid_Dierckx[x] = evaluate(V_1D_int_Dierckx, xt)
  V_1D_actual[x] = f(xt)
end

for y = 1:5*ypoints
  for z = 1:5*zpoints
    yt = yoffgrid[y]
    zt = zoffgrid[z]
    V_2D_offgrid_Lin[y, z] = V_2D_int_Lin[yt, zt]
    V_2D_offgrid_Quad[y, z] = V_2D_int_Quad[yt, zt]
    V_2D_offgrid_Dierckx[y, z] = evaluate(V_2D_int_Dierckx, yt, zt)
    V_2D_actual[y, z] = f(yt + zt)
  end
end

######################################################################

# Errors
Err_1D_Lin = V_1D_actual - V_1D_offgrid_Lin
Err_1D_Quad = V_1D_actual - V_1D_offgrid_Quad
Err_1D_Dierckx = V_1D_actual - V_1D_offgrid_Dierckx
Err_2D_Lin = V_2D_actual - V_2D_offgrid_Lin
Err_2D_Quad = V_2D_actual - V_2D_offgrid_Quad
Err_2D_Dierckx = V_2D_actual - V_2D_offgrid_Dierckx

sum(abs(Err_1D_Lin))
sum(abs(Err_1D_Quad))
sum(abs(Err_1D_Dierckx))

sum(abs(Err_2D_Lin))
sum(abs(Err_2D_Quad))
sum(abs(Err_2D_Dierckx))

######################################################################

# Plots
fig, ax = PyPlot.subplots(3,1)
ax[1,1][:plot](V_1D_actual[1:20], label = "True Value")
ax[1,1][:plot](V_1D_offgrid_Lin[1:20], label = "Linearly Interpolated Value")
ax[1,1][:plot](Err_1D_Lin[1:20], label = "Grid Linear Interpolation Error")
ax[1,1][:legend]()
ax[2,1][:plot](V_1D_actual[1:20], label = "True Value")
ax[2,1][:plot](V_1D_offgrid_Quad[1:20], label = "Quadratically Interpolated Value")
ax[2,1][:plot](Err_1D_Quad[1:20], label = "Grid Quadratic Interpolation Error")
ax[2,1][:legend]()
ax[3,1][:plot](V_1D_actual[1:20], label = "True Value")
ax[3,1][:plot](V_1D_offgrid_Dierckx[1:20], label = "Quadratically Interpolated Value")
ax[3,1][:plot](Err_1D_Dierckx[1:20], label = "Dierckx Interpolation Error")
ax[3,1][:legend]()
fig[:suptitle]("Interpolation in 1D")
plt.show()

fig = figure()
ax = fig[:add_subplot](221, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_actual[1:ypoints, 1:zpoints, 1]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Actual Value")

ax = fig[:add_subplot](222, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Lin[1:ypoints, 1:zpoints, 1]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Grid Linear Interpolation")

ax = fig[:add_subplot](223, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Quad[1:ypoints, 1:zpoints, 1]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Grid Quadratic Interpolation")

ax = fig[:add_subplot](224, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Dierckx[1:20, 1:15, 1]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Dierckx Spline Interpolation")
plt.show()
