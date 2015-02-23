################################################################################

using  ApproXD, Dierckx, Grid, NumericalMath, PyPlot, PyCall
@pyimport seaborn as sns

################################################################################
# Interpolation schemes to consider:
methods_1D = ["Grid (linear)", "Grid (quadratic)", "ApproXD (linear)",
              "ApproXD (linear, irregular grid)", "Dierckx (spline)",
              "Dierckx (spline, irregular grid)", "NumericalMath(pchip)"]

methods_2D = ["Grid (linear)", "Grid (quadratic)", "ApproXD (linear)",
              "ApproXD (linear, irregular grid)", "Dierckx (spline)",
              "Dierckx (spline, irregular grid)"]

# Function to interpolate
f(x::Float64) = -(x^(-1.0))

# Grids for interpolation (regular/irregular)
zmin = -0.1
zmax = 10.0
ymin = 0.2
ymax = 30.0
xpoints = 100
ypoints = 20
zpoints = 15

xmin = ymin + zmin
xmax = ymax + zmax

xgrid = linspace(xmin, xmax, xpoints)
ygrid = linspace(ymin, ymax, ypoints)
zgrid = linspace(zmin, zmax, zpoints)

power = 2.0  # Curvature of the irregular grids
x_irrgrid = Array(Float64, xpoints)
x_irrgridexp = Array(Float64, xpoints)
xdistexp = (xmax - xmin)^(1/power)
xinc = xdistexp/(xpoints-1)
for i = 1:xpoints
  x_irrgridexp[i] = (i-1)*xinc
end
x_irrgrid[:] = x_irrgridexp.^power + xmin

y_irrgrid = Array(Float64, ypoints)
y_irrgridexp = Array(Float64, ypoints)
ydisteyp = (ymax - ymin)^(1/power)
yinc = ydisteyp/(ypoints-1)
for i = 1: ypoints
  y_irrgridexp[i] = (i-1)*yinc
end
y_irrgrid[:] = y_irrgridexp.^power + ymin

z_irrgrid = Array(Float64, zpoints)
z_irrgridexp = Array(Float64, zpoints)
zdistezp = (zmax - zmin)^(1/power)
zinc = zdistezp/(zpoints-1)
for i = 1: zpoints
  z_irrgridexp[i] = (i-1)*zinc
end
z_irrgrid[:] = z_irrgridexp.^power + zmin

################################################################################
# Evaluate functions on grid points

V_1D = Array(Float64, xpoints)
V_1D_irr = similar(V_1D)
V_2D = Array(Float64, (ypoints, zpoints))
V_2D_irr = similar(V_2D)

# Values on grid, 1 dimension
for i = 1:xpoints
  V_1D[i] = f(xgrid[i])
  V_1D_irr[i] = f(x_irrgrid[i])
end

# Values on grid, 2 dimensions
for i = 1:ypoints
  for j = 1:zpoints
    V_2D[i, j] = f(ygrid[i] + zgrid[j])
    V_2D_irr[i, j] = f(y_irrgrid[i] + z_irrgrid[j])
  end
end

################################################################################
# Interpolate in 1 dimension

# Grid (regular grid only)
xrange = range(xgrid[1], xgrid[2]-xgrid[1], xpoints)

V_1D_int_Lin =
  CoordInterpGrid(xrange, V_1D, BCnearest, InterpLinear)
V_1D_int_Quad =
  CoordInterpGrid(xrange, V_1D, BCnearest, InterpQuadratic)

# ApproXD (regular/irregular grid)
xs = Array{Float64, 1}[]
push!(xs, xgrid)
V_1D_int_ApproXD = lininterp(V_1D, xs)
xsirr = Array{Float64, 1}[]
push!(xsirr, x_irrgrid)
V_1D_int_irr_ApproXD = lininterp(V_1D_irr, xsirr)

# Dierckx (regular/irregular grid)
V_1D_int_Dierckx = Spline1D(xgrid, V_1D)
V_1D_int_irr_Dierckx = Spline1D(x_irrgrid, V_1D_irr)

# NumericalMath - Piecewise cubic hermitean interpolating Polynomial
# As in Matlab, NumericalMath doesn't separate interpolation and
# evaluation, there is only one call to pchip(xgrid, V_1D, x)

################################################################################
# Interpolate in 2 dimensions

# Grid
yrange = range(ygrid[1], ygrid[2]-ygrid[1], ypoints)
zrange = range(zgrid[1], zgrid[2]-zgrid[1], zpoints)

V_2D_int_Lin =
  CoordInterpGrid((yrange, zrange), V_2D, BCnearest, InterpLinear)
V_2D_int_Quad =
  CoordInterpGrid((yrange, zrange), V_2D, BCnearest, InterpQuadratic)
#V_2D_int_Cub = CoordInterpGrid((yrange, zrange), V_2D, BCnan, InterpCubic)

# ApproXD (regular/irregular grid)
yz = Array{Float64, 1}[]
push!(yz, ygrid)
push!(yz, zgrid)
V_2D_int_ApproXD = lininterp(V_2D, yz)
yz_irr = Array{Float64, 1}[]
push!(yz_irr, y_irrgrid)
push!(yz_irr, z_irrgrid)
V_2D_int_irr_ApproXD = lininterp(V_2D_irr, yz_irr)

# Dierckx (regular/irregular grid)
V_2D_int_Dierckx = Spline2D(ygrid, zgrid, V_2D)
V_2D_int_irr_Dierckx = Spline2D(y_irrgrid, z_irrgrid, V_2D_irr)

################################################################################
# Evaluate Regular and Irregular Grid Interpolants Off-Grid

xop = 20*xpoints  # Number of off-grid points
yop = 20*ypoints
zop = 20*zpoints

xoffgrid = linspace(xgrid[1], xgrid[end], xop)
yoffgrid = linspace(ygrid[1], ygrid[end], yop)
zoffgrid = linspace(zgrid[1], zgrid[end], zop)

V_1D_offgrid = Array(Float64, (length(methods_1D), xop))
V_1D_actual = similar(V_1D_offgrid)

V_2D_offgrid = Array(Float64, (length(methods_2D), yop, zop))
V_2D_actual = similar(V_2D_offgrid)

for i = 1:xop
  # Grid
  V_1D_offgrid[1, i] = V_1D_int_Lin[xoffgrid[i]]
  V_1D_offgrid[2, i] = V_1D_int_Quad[xoffgrid[i]]
  # ApproXD (regular/irregular grid)
  V_1D_offgrid[3, i] = getValue(V_1D_int_ApproXD, [xoffgrid[i]])[1]
  V_1D_offgrid[4, i] = getValue(V_1D_int_irr_ApproXD, [xoffgrid[i]])[1]
  # Dierckx (regular/irregular grid)
  V_1D_offgrid[5, i] = evaluate(V_1D_int_Dierckx, xoffgrid[i])
  V_1D_offgrid[6, i] = evaluate(V_1D_int_irr_Dierckx, xoffgrid[i])
  # NumericalMath
  V_1D_offgrid[7, i] = pchip(xgrid, V_1D, xoffgrid[i])
  # Actual Function value
  V_1D_actual[:, i] = f(xoffgrid[i])
end

for i = 1:yop
  for j = 1:zop
    yt = yoffgrid[i]
    zt = zoffgrid[j]
    # Grid
    V_2D_offgrid[1, i, j] = V_2D_int_Lin[yt, zt]
    V_2D_offgrid[2, i, j] = V_2D_int_Quad[yt, zt]
    # ApproXD (regular/irregular grid)
    V_2D_offgrid[3, i, j] = getValue(V_2D_int_ApproXD, [yt, zt])[1]
    V_2D_offgrid[4, i, j] = getValue(V_2D_int_irr_ApproXD, [yt, zt])[1]
    # Dierckx (regular/irregular grid)
    V_2D_offgrid[5, i, j] = evaluate(V_2D_int_Dierckx, yt, zt)
    V_2D_offgrid[6, i, j] = evaluate(V_2D_int_irr_Dierckx, yt, zt)
    # Actual Function value
    V_2D_actual[:, i, j] = f(yt + zt)
  end
end

######################################################################

# Calculate Errors
Error_1D = (V_1D_actual - V_1D_offgrid)./V_1D_actual
Error_2D = (V_2D_actual - V_2D_offgrid)./V_2D_actual

######################################################################
# Print results
@printf "\n∑|f(x)-f(x_int)/f(x)|/n in 1 dimension\n"
for i = 1:length(methods_1D)
  @printf "\t%s: %.4f\n" methods_1D[i] sum(abs(Error_1D[i, :]))/xop
end

@printf "\n∑|f(x)-f(x_int)/f(x)|/n in 2 dimensions\n"
for i = 1:length(methods_2D)
  @printf "\t%s: %.4f\n" methods_2D[i] sum(abs(Error_2D[i, :, :]))/(yop*zop)
end

# Plot 1D results
fig, ax = PyPlot.subplots(2,1, figsize = (12,10))
for i = 1:7
  ax[1,1][:plot](V_1D_offgrid[i, 1:end/20]', label = methods_1D[i])
  ax[2,1][:plot](Error_1D[i, 1:end/20]')
end
ax[1,1][:plot](V_1D_actual[1, 1:end/20]', label = "True Value",
               linestyle = ":", color = "black")
ax[1,1][:set_title]("True Value vs. Interpolated Values")
ax[1,1][:legend](loc = "best")
ax[2,1][:set_title]("Relative Error")
ax[2,1][:legend](loc = "best")
fig[:suptitle]("Interpolation in One Dimension", fontsize = 16)
