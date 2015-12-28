################################################################################
#################### Interpolations for value and policy functions #############
################################################################################

using Interpolations

################################################################################

function interpolateV{T<:AbstractFloat}(v::Array{T,5}, wgrid::Array{T,1},
  hgrid::Array{T,1}, agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1})

  interpolate((wgrid, hgrid, agrid, bgrid, zgrid), v, Gridded(Linear()))
end

################################################################################

function interpolateV{T<:AbstractFloat}(v::Array{T,3}, xgrid::Array{T,1},
                      ygrid::Array{T,1}, hgrid::Array{T,1})

  interpolate((xgrid, ygrid, hgrid), v, Gridded(Linear()))
end

################################################################################
