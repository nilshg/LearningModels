################################################################################
#################### Interpolations for value and policy functions #############
################################################################################

using Interpolations

################################################################################

function interpolateV{T<:AbstractFloat}(v::Array{T,4}, xgrid::Array{T,1},
            agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1})

  return interpolate((xgrid, agrid, bgrid, zgrid), v, Gridded(Linear()))
end

################################################################################

function interpolateV(v::Array{Float64,2}, xgrid::Array{Float64,1},
                      ygrid::Array{Float64,1})

  return interpolate((xgrid, ygrid), v, Gridded(Linear()))
end

################################################################################

function get_lins{T<:AbstractFloat}(v::Array{T,4}, xgrid::Array{T,1},
  agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1})

  lin_int = Array(Lininterp, length(xgrid))
  for x = 1:length(xgrid)
    lin_int[x] = interpolateV(squeeze(v[x,:,:,:], 1), agrid, bgrid, zgrid)
  end
  return lin_int
end
