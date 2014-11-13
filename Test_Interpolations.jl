######################################################################
using PyPlot
using Grid
using Dierckx
using NumericalMath
#using ApproXD

meshgrid(x::Vector, y::Vector) = (repmat(x, 1, length(y))',
                                  repmat(y, 1, length(x)))

######################################################################

# Function to interpolate 
f(x::Float64) = -(x^(-1.0))

# Grids for interpolation (regular/irregular)
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

power = 3  # Curvature of the irregular grids
x_irrgrid = zeros(xpoints)
x_irrgridexp = zeros(xpoints)
xdistexp = (xmax - xmin)^(1/power)
xinc = xdistexp/(xpoints-1)
for i = 1: xpoints
  x_irrgridexp[i] = (i-1)*xinc
end
x_irrgrid[:] = x_irrgridexp.^power + xmin

y_irrgrid = zeros(ypoints)
y_irrgridexp = zeros(ypoints)
ydisteyp = (ymax - ymin)^(1/power)
yinc = ydisteyp/(ypoints-1)
for i = 1: ypoints
  y_irrgridexp[i] = (i-1)*yinc
end
y_irrgrid[:] = y_irrgridexp.^power + ymin

z_irrgrid = zeros(zpoints)
z_irrgridexp = zeros(zpoints)
zdistezp = (zmax - zmin)^(1/power)
zinc = zdistezp/(zpoints-1)
for i = 1: zpoints
  z_irrgridexp[i] = (i-1)*zinc
end
z_irrgrid[:] = z_irrgridexp.^power + zmin

######################################################################
# Evaluate functions on grid points

V_1D = zeros(xpoints)
V_1D_irr = zeros(xpoints)
V_2D = zeros(ypoints, zpoints)
V_2D_irr = zeros(ypoints, zpoints)

# Values on grid, 1 dimension
for x = 1:xpoints
  xt = xgrid[x]
  xt2 = x_irrgrid[x]
  V_1D[x] = f(xt)
  V_1D_irr[x] = f(xt2)
end

# Values on grid, 2 dimensions
for z = 1:zpoints
  for y = 1:ypoints
    zt = zgrid[z]
    yt = ygrid[y]
    zt2 = z_irrgrid[z]
    yt2 = y_irrgrid[y]
    V_2D[y, z, 1] = f(yt+zt)
    V_2D_irr[y, z, 1] = f(yt2 + zt2)
  end
end

######################################################################
# Interpolate in 1 dimension

xrange = range(xgrid[1], xgrid[2]-xgrid[1], xpoints)

# Grid
V_1D_int_Lin = CoordInterpGrid(xrange, V_1D, BCnearest, InterpLinear)
V_1D_int_Quad = CoordInterpGrid(xrange, V_1D, BCnearest, InterpQuadratic)
#V_1D_int_Cub = CoordInterpGrid(xrange, V_1D[:, 1], BCnil, InterpCubic)

# ApproXD
# V_1D_int_ApproXD = 

# Dierckx (regular/irregular grid)
V_1D_int_Dierckx = Spline1D(xgrid, V_1D)
V_1D_int_irr_Dierckx = Spline1D(x_irrgrid, V_1D_irr)

######################################################################  
# Interpolate in 2 dimensions

yrange = range(ygrid[1], ygrid[2]-ygrid[1], ypoints)
zrange = range(zgrid[1], zgrid[2]-zgrid[1], zpoints)

# Grid
V_2D_int_Lin = CoordInterpGrid((yrange, zrange), V_2D, BCnearest, InterpLinear)
V_2D_int_Quad = CoordInterpGrid((yrange, zrange), V_2D, BCnearest, InterpQuadratic)
#V_2D_int_Cub = CoordInterpGrid((yrange, zrange), V_2D, BCnan, InterpCubic)

# ApproXD
# V_2D_int_ApproXD = 

# Dierckx (regular/irregular grid)
V_2D_int_Dierckx = Spline2D(ygrid, zgrid, V_2D)
V_2D_int_irr_Dierckx = Spline2D(y_irrgrid, z_irrgrid, V_2D_irr)

# NumericalMath - Piecewise cubic hermitean interpolating Polynomial
# As in Matlab, NumericalMath doesn't separate interpolation and 
# evaluation, there is only one call to pchip(xgrid, V_1D, x)

######################################################################
# Evaluate Regular and Irregular Grid Interpolants Off-Grid

xop = 5*xpoints  # Number of off-grid points 
yop = 5*ypoints
zop = 5*zpoints

xoffgrid = linspace(xgrid[1], xgrid[end], xop)
yoffgrid = linspace(ygrid[1], ygrid[end], yop)
zoffgrid = linspace(zgrid[1], zgrid[end], zop)

V_1D_offgrid_Lin = zeros(xop)
V_1D_offgrid_Quad = zeros(xop)
V_1D_offgrid_Dierckx = zeros(xop)
V_1D_offgrid_pchip = zeros(xop)
V_1D_irr_offgrid_Dierckx = zeros(xop)
V_1D_actual = zeros(xop)

V_2D_offgrid_Lin =  zeros(yop, zop)
V_2D_offgrid_Quad = zeros(yop, zop)
V_2D_offgrid_Dierckx = zeros(yop, zop)
V_2D_irr_offgrid_Dierckx = zeros(yop, zop)
V_2D_actual = zeros(yop, zop)

for x = 1:xop
  xt = xoffgrid[x]
  # Grid
  V_1D_offgrid_Lin[x] = V_1D_int_Lin[xt]
  V_1D_offgrid_Quad[x] = V_1D_int_Quad[xt]
  # Dierckx (regular/irregular grid)
  V_1D_offgrid_Dierckx[x] = evaluate(V_1D_int_Dierckx, xt)
  V_1D_irr_offgrid_Dierckx[x] = evaluate(V_1D_int_irr_Dierckx, xt)
  # NumericalMath
  V_1D_offgrid_pchip[x] = pchip(xgrid, V_1D, xt)
  # Actual Function value
  V_1D_actual[x] = f(xt)
end

for y = 1:yop
  for z = 1:zop
    yt = yoffgrid[y]
    zt = zoffgrid[z]
    # Grid
    V_2D_offgrid_Lin[y, z] = V_2D_int_Lin[yt, zt]
    V_2D_offgrid_Quad[y, z] = V_2D_int_Quad[yt, zt]
    # Dierckx (regular/irregular grid)
    V_2D_offgrid_Dierckx[y, z] = evaluate(V_2D_int_Dierckx, yt, zt)
    V_2D_irr_offgrid_Dierckx[y, z] = evaluate(V_2D_int_irr_Dierckx, yt, zt)
    # Actual Function value 
    V_2D_actual[y, z] = f(yt + zt)
  end
end

######################################################################
# Calculate Errors

Err_1D_Lin = (V_1D_actual - V_1D_offgrid_Lin)./V_1D_actual
Err_1D_Quad = (V_1D_actual - V_1D_offgrid_Quad)./V_1D_actual
Err_1D_Dierckx = (V_1D_actual - V_1D_offgrid_Dierckx)./V_1D_actual
Err_1D_pchip = (V_1D_actual - V_1D_offgrid_pchip)./V_1D_actual
Err_1D_irr_Dierckx = (V_1D_actual - V_1D_irr_offgrid_Dierckx)./V_1D_actual
Err_2D_Lin = (V_2D_actual - V_2D_offgrid_Lin)./V_2D_actual
Err_2D_Quad = (V_2D_actual - V_2D_offgrid_Quad)./V_2D_actual
Err_2D_Dierckx = (V_2D_actual - V_2D_offgrid_Dierckx)./V_2D_actual
Err_2D_irr_Dierckx = (V_2D_actual - V_2D_irr_offgrid_Dierckx)./V_2D_actual

@printf "∑(abs(ɛ)) for 1D linear interpolation with Grid is %.3f\n" sum(abs(Err_1D_Lin))/xop
@printf "∑(abs(ɛ)) for 1D quadratic interpolation with Grid is %.3f\n" sum(abs(Err_1D_Quad))/xop
@printf "∑(abs(ɛ)) for 1D spline interpolation on a regular grid with Dierckx is %.3f\n" sum(abs(Err_1D_Dierckx))/xop
@printf "∑(abs(ɛ)) for 1D spline interpolation on an irregular grid with Dierckx is %.3f\n" sum(abs(Err_1D_irr_Dierckx))/xop
@printf "∑(abs(ɛ)) for 1D pchip interpolation with NumericalMath is %.2f\n\n" sum(abs(Err_1D_pchip))/xop

@printf "∑(abs(ɛ)) for 2D linear interpolation with Grid is %.3f\n" sum(abs(Err_2D_Lin))/(yop*zop)
@printf "∑(abs(ɛ)) for 2D quadratic interpolation with Grid is %.3f\n" sum(abs(Err_2D_Quad))/(yop*zop)
@printf "∑(abs(ɛ)) for 2D spline interpolation on a regular grid with Dierckx is %.3f\n" sum(abs(Err_2D_Dierckx))/(yop*zop)
@printf "∑(abs(ɛ)) for 2D spline interpolation on an irregular grid with Dierckx is %.3f\n" sum(abs(Err_2D_irr_Dierckx))/(yop*zop)

######################################################################

# Plots
fig, ax = PyPlot.subplots(2,3)
ax[1,1][:plot](V_1D_actual[1:20], label = "True Value")
ax[1,1][:plot](V_1D_offgrid_Lin[1:20], label = "Grid Linearly Interpolated Value")
ax[1,1][:plot](Err_1D_Lin[1:20], label = "Interpolation Error")
ax[1,1][:legend](loc = 0)
ax[2,1][:plot](V_1D_actual[1:20], label = "True Value")
ax[2,1][:plot](V_1D_offgrid_Quad[1:20], label = "Grid Quadratically Interpolated Value")
ax[2,1][:plot](Err_1D_Quad[1:20], label = "Interpolation Error")
ax[2,1][:legend](loc = 0)
ax[1,2][:plot](V_1D_actual[1:20], label = "True Value")
ax[1,2][:plot](V_1D_offgrid_Dierckx[1:20], label = "Dierckx Interpolated Value (regular grid)")
ax[1,2][:plot](Err_1D_Dierckx[1:20], label = "Interpolation Error")
ax[1,2][:legend](loc = 0)
ax[2,2][:plot](V_1D_actual[1:20], label = "True Value")
ax[2,2][:plot](V_1D_offgrid_pchip[1:20], label = "PCHIP Interpolated Value")
ax[2,2][:plot](Err_1D_pchip[1:20], label = "Interpolation Error")
ax[2,2][:legend](loc = 0)
ax[2,3][:plot](V_1D_actual[1:20], label = "True Value")
ax[2,3][:plot](V_1D_irr_offgrid_Dierckx[1:20], label = "Dierckx Interpolated Value (irregular grid)")
ax[2,3][:plot](Err_1D_irr_Dierckx[1:20], label = "Interpolation Error")
ax[2,3][:legend](loc = 0)
fig[:suptitle]("Interpolation in 1D")
plt.show()

fig2 = figure()
ax = fig2[:add_subplot](221, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Lin[1:ypoints, 1:zpoints]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Grid Linear Interpolation")

ax = fig2[:add_subplot](222, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Quad[1:ypoints, 1:zpoints]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Grid Quadratic Interpolation")

ax = fig2[:add_subplot](223, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_offgrid_Dierckx[1:20, 1:15]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Dierckx Spline Interpolation (regular grid")

ax = fig2[:add_subplot](224, projection="3d")
yg, zg = meshgrid(yoffgrid[1:ypoints], zoffgrid[1:zpoints])
vg = V_2D_irr_offgrid_Dierckx[1:ypoints, 1:zpoints]'
ax[:plot_surface](yg, zg, vg, rstride = 1, cstride = 1,
                    cmap=ColorMap("jet"), alpha=0.5, linewidth=0.25)
ax[:set_title]("Dierckx Spline Interpolation (irregular grid)")
plt.show()
