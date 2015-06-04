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

function interpolateV(v::Array{Float64, 4}, xgrid::Array{Float64, 1},
                        agrid::Array{Float64, 1}, bgrid::Array{Float64, 1},
                        zgrid::Array{Float64, 1})

  xabz = Array{Float64, 1}[]
  push!(xabz, xgrid)
  push!(xabz, agrid)
  push!(xabz, bgrid)
  push!(xabz, zgrid)
  return Lininterp(v, xabz)
end

#######################################################################################

function interpolateV(v::Array{Float64, 2}, xgrid::Array{Float64, 1},
                        ygrid::Array{Float64, 1})

  xy = Array{Float64, 1}[]
  push!(xy, xgrid)
  push!(xy, ygrid)
  return Lininterp(v, xy)
end
