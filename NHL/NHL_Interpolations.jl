################################################################################
#################### Interpolations for value and policy functions #############
################################################################################

using ApproXD

################################################################################

function interpolateV{T<:AbstractFloat}(v::Array{T,4}, xgrid::Array{T,1},
            agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1})

  xabz = Array{Float64,1}[]
  push!(xabz, xgrid); push!(xabz, agrid); push!(xabz, bgrid); push!(xabz, zgrid)
  return Lininterp(v, xabz)
end

################################################################################

function interpolateV(v::Array{Float64,2}, xgrid::Array{Float64,1},
                      ygrid::Array{Float64,1})

  xy = Array{Float64,1}[]
  push!(xy, xgrid); push!(xy, ygrid)
  return Lininterp(v, xy)
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
