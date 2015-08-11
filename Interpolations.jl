#######################################################################################
################## Interpolations - Interpolation for value Functions #################
#######################################################################################

# Contains:
# interpolateV(V::Array{Float64, 5}, wgrid, ygrid, agrid, bgrid, zgrid, t):
#     Interpolate value function w/o habits using CoordInterpGrid
# interpolateV(V::Array{Float64, 3}, wgrid, hgrid, ygrid, t):
#     Interpolate retirement value w/o habits function using CoordInterpGrid function

#######################################################################################

using ApproXD

#######################################################################################

function interpolateV(v::Array{Float64,4}, xgrid::Array{Float64,1},
                      agrid::Array{Float64,1}, bgrid::Array{Float64,1},
                      zgrid::Array{Float64,1})

  xabz = Array{Float64,1}[]
  push!(xabz, xgrid)
  push!(xabz, agrid)
  push!(xabz, bgrid)
  push!(xabz, zgrid)
  return Lininterp(v, xabz)
end

#######################################################################################

function interpolateV(v::Array{Float64,2}, xgrid::Array{Float64,1},
                      ygrid::Array{Float64,1})

  xy = Array{Float64,1}[]
  push!(xy, xgrid)
  push!(xy, ygrid)
  return Lininterp(v, xy)
end

#######################################################################################

function interpolateV(v::Array{Float64, 5}, wgrid::Array{Float64,1},
                      hgrid::Array{Float64,1}, agrid::Array{Float64,1},
                      bgrid::Array{Float64,1}, zgrid::Array{Float64,1})

  wr = range(wgrid[1], wgrid[2] - wgrid[1], size(wgrid, 1))
  hr = range(hgrid[1], hgrid[2] - hgrid[1], size(hgrid, 1))
  ar = range(agrid[1], agrid[2] - agrid[1], size(agrid, 1))
  br = range(bgrid[1], bgrid[2] - bgrid[1], size(bgrid, 1))
  zr = range(zgrid[1], zgrid[2] - zgrid[1], size(zgrid, 1))

  CoordInterpGrid((wr, hr, ar, br, zr), v, BCnearest, InterpLinear)
end

#######################################################################################

function interpolateV(v::Array{Float64,3}, xgrid::Array{Float64,1},
                      ygrid::Array{Float64,1}, hgrid::Array{Float64,1})

  xyh = Array{Float64,1}[]
  push!(xyh, xgrid)
  push!(xyh, ygrid)
  push!(xyh, hgrid)
  return Lininterp(v, xyh)
end

#######################################################################################
