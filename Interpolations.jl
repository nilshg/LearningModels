#######################################################################################
################## Interpolations - Interpolation for value Functions #################
#######################################################################################

# Contains:
# interpolatev(v::Array{Float64, 6}, wgrid, ygrid, hgrid, agrid, bgrid, zgrid, t):
#     Interpolate value function with habits using Grid.jl's CoordInterpGrid function
# interpolatev(v::Array{Float64, 4}, wgrid, hgrid, ygrid, t):
#     Interpolate retirement value function with habits using CoordInterpGrid
# interpolateV(V::Array{Float64, 5}, wgrid, ygrid, agrid, bgrid, zgrid, t):
#     Interpolate value function w/o habits using CoordInterpGrid
# interpolateV(V::Array{Float64, 3}, wgrid, hgrid, ygrid, t):
#     Interpolate retirement value w/o habits function using CoordInterpGrid function

#######################################################################################

using Grid

#######################################################################################

function interpolatev(v::Array{Float64, 6},     wgrid::Array{Float64, 2},
                      hgrid::Array{Float64, 2}, agrid::Array{Float64, 1},
                      bgrid::Array{Float64, 1}, zgrid::Array{Float64, 1}, t::Int64)

  wrange = range(wgrid[1, t], wgrid[2, t] - wgrid[1, t], size(wgrid, 1))
  hrange = range(hgrid[1, t], hgrid[2, t] - hgrid[1, t], size(hgrid, 1))
  arange = range(agrid[1], agrid[2] - agrid[1], size(agrid, 1))
  brange = range(bgrid[1], bgrid[2] - bgrid[1], size(bgrid, 1))
  zrange = range(zgrid[1], zgrid[2] - zgrid[1], size(zgrid, 1))

  CoordInterpGrid((wrange, hrange, arange, brange, zrange), v[:, :, :, :, :, t],
                  BCnearest, InterpLinear)
end

#######################################################################################

function interpolatev(v::Array{Float64, 4}, wgrid::Array{Float64, 2},
                      hgrid::Array{Float64, 2}, ygrid::Array{Float64, 1}, t::Int64)

  wrange = range(wgrid[1, t], wgrid[2, t] - wgrid[1, t], size(wgrid, 1))
  hrange = range(hgrid[1, t], hgrid[2, t] - hgrid[1, t], size(hgrid, 1))
  yrange = range(ygrid[1], ygrid[2] - ygrid[1], size(ygrid, 1))

  CoordInterpGrid((wrange, hrange, yrange), v[:, :, :, t], BCnearest, InterpLinear)
end

#######################################################################################

function interpolatev(v::Array{Float64, 5}, wgrid::Array{Float64, 2},
                      agrid::Array{Float64, 1}, bgrid::Array{Float64, 1},
                      zgrid::Array{Float64, 1}, t::Int64)

  wrange = range(wgrid[1, t], wgrid[2, t] - wgrid[1, t], size(wgrid, 1))
  arange = range(agrid[1], agrid[2] - agrid[1], size(agrid, 1))
  brange = range(bgrid[1], bgrid[2] - bgrid[1], size(bgrid, 1))
  zrange = range(zgrid[1], zgrid[2] - zgrid[1], size(zgrid, 1))

  CoordInterpGrid((wrange, arange, brange, zrange), v[:, :, :, :, t], BCnearest,
                  InterpLinear)
end

#######################################################################################

function interpolatev(v::Array{Float64,3}, wgrid::Array{Float64,2},
                        ygrid::Array{Float64,1}, t::Int64)

  wrange = range(wgrid[1, t], wgrid[2, t] - wgrid[1, t], size(wgrid, 1))
  yrange = range(ygrid[1], ygrid[2]-ygrid[1], size(ygrid, 1))

  CoordInterpGrid((wrange, yrange), v[:, :, t], BCnearest, InterpLinear)
end

#######################################################################################

function interpolatev(v::Array{Float64, 2}, xgrid::Array{Float64, 2}, t::Int64)

  xrange = range(xgrid[1, t], xgrid[2, t] - xgrid[1, t], size(xgrid, 1))

  CoordInterpGrid(xrange, v[:, t], BCnearest, InterpLinear)
end
