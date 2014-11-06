#######################################################################################
################## Interpolations - Interpolation for Value Functions #################
#######################################################################################

# Contains:
# interpolateV(V::Array{Float64, 6}, wgrid, ygrid, hgrid, agrid, bgrid, zgrid, t):
#     Interpolate value function with habits using Grid.jl's CoordInterpGrid function
# interpolateV(V::Array{Float64, 4}, wgrid, hgrid, ygrid, t):
#     Interpolate retirement value function with habits using CoordInterpGrid
# interpolateV(V::Array{Float64, 5}, wgrid, ygrid, agrid, bgrid, zgrid, t):
#     Interpolate value function w/o habits using CoordInterpGrid
# interpolateV(V::Array{Float64, 3}, wgrid, hgrid, ygrid, t):
#     Interpolate retirement value w/o habits function using CoordInterpGrid function

#######################################################################################

function interpolateV(V::Array{Float64, 6}, wgrid::Array, hgrid::Array,
                      agrid::Array, bgrid::Array, zgrid::Array, t::Int64)

  wrange = wgrid[end, t] - wgrid[1, t]
  wstep = wrange / (size(wgrid, 1)-1)
  wrange = wgrid[1, t]:wstep:wgrid[end, t]

  if length(wrange) < size(wgrid, 1)
    wrange = wgrid[1, t]:wstep:(wgrid[end, t]+0.0001)
  end

  hrange = hgrid[end, t] - hgrid[1, t]
  hstep = hrange / (size(hgrid, 1)-1)
  hrange = hgrid[1, t]:hstep:hgrid[end, t]

  if length(hrange) < size(hgrid, 1)
    hrange = hgrid[1, t]:hstep:(hgrid[end, t]+0.0001)
  end

  arange = agrid[end, t] - agrid[1, t]
  astep = arange / (size(agrid, 1)-1)
  arange = agrid[1, t]:astep:agrid[end, t]

  if length(arange) < size(agrid, 1)
    arange = agrid[1, t]:astep:(agrid[end, t]+0.0001)
  end

  brange = bgrid[end, t] - bgrid[1, t]
  bstep = brange / (size(bgrid, 1)-1)
  brange = bgrid[1, t]:bstep:bgrid[end, t]

  if length(brange) < size(bgrid, 1)
    brange = bgrid[1, t]:bstep:(bgrid[end, t]+0.0001)
  end

  zrange = zgrid[end, t] - zgrid[1, t]
  zstep = zrange / (size(zgrid, 1)-1)
  zrange = zgrid[1, t]:zstep:zgrid[end, t]

  if length(zrange) < size(zgrid, 1)
    zrange = zgrid[1, t]:zstep:(zgrid[end, t]+0.0001)
  end

  CoordInterpGrid((wrange, hrange, arange, brange, zrange), V[:, :, :, :, :, t],
                  BCnearest, InterpQuadratic)
end

#######################################################################################

function interpolateV(V::Array{Float64, 4}, wgrid_R::Array, hgrid_R::Array, 
                      ygrid_R::Array, t::Int64)

  wrange = wgrid_R[end, t] - wgrid_R[1, t]
  wstep = wrange / (size(wgrid_R, 1)-1)
  wrange = wgrid_R[1, t]:wstep:wgrid_R[end, t]

  if length(wrange) < size(wgrid_R, 1)
    wrange = wgrid_R[1, t]:wstep:(wgrid_R[end, t]+0.0001)
  end

  hrange = hgrid_R[end, t] - hgrid_R[1, t]
  hstep = hrange / (size(hgrid_R, 1)-1)
  hrange = hgrid_R[1, t]:hstep:hgrid_R[end, t]

  yrange = ygrid_R[end] - ygrid_R[1]
  ystep = yrange / (size(ygrid_R, 1)-1)
  yrange = ygrid_R[1]:ystep:ygrid_R[end]

  CoordInterpGrid((wrange, hrange, yrange), V[:, :, :, t], BCnearest, InterpLinear)

end

#######################################################################################

function interpolateV(V::Array{Float64, 5}, wgrid::Array{Float64, 2},
                      agrid::Array{Float64, 2}, bgrid::Array{Float64, 2},
                      zgrid::Array{Float64, 2}, t::Int64)

  wrange = wgrid[end, t] - wgrid[1, t]
  wstep = wrange / (size(wgrid, 1)-1)
  wrange = wgrid[1, t]:wstep:wgrid[end, t]

  if length(wrange) < size(wgrid, 1)
    wrange = wgrid[1, t]:wstep:(wgrid[end, t] + 0.0001)
  end

  arange = agrid[end, t] - agrid[1, t]
  astep = arange / (size(agrid, 1)-1)
  arange = agrid[1, t]:astep:agrid[end, t]

  if length(arange) < size(agrid, 1)
    arange = agrid[1, t]:astep:(agrid[end, t] + 0.0001)
  end

  brange = bgrid[end, t] - bgrid[1, t]
  bstep = brange / (size(bgrid, 1)-1)
  brange = bgrid[1, t]:bstep:bgrid[end, t]

  if length(brange) < size(bgrid, 1)
    brange = bgrid[1, t]:bstep:(bgrid[end, t] + 0.0001)
  end

  zrange = zgrid[end, t] - zgrid[1, t]
  zstep = zrange / (size(zgrid, 1)-1)
  zrange = zgrid[1, t]:zstep:zgrid[end, t]

  if length(zrange) < size(zgrid, 1)
    zrange = zgrid[1, t]:zstep:(zgrid[end, t] + 0.0001)
  end

  CoordInterpGrid((wrange, arange, brange, zrange), V[:, :, :, :, t], BCnearest,
                  InterpQuadratic)
end

#######################################################################################

function interpolateV(V::Array{Float64,3}, wgrid::Array{Float64,2},
                        ygrid::Array{Float64,1}, t::Int64)

  wrange = wgrid[end, t] - wgrid[1, t]
  wstep = wrange / (size(wgrid, 1)-1)
  wrange = wgrid[1, t]:wstep:wgrid[end, t]

  if length(wrange) < size(wgrid, 1)
    wrange = wgrid[1, t]:wstep:(wgrid[end, t] + 0.0001)
  end

  yrange = ygrid[end] - ygrid[1]
  ystep = yrange / (size(ygrid, 1)-1)
  yrange = ygrid[1]:ystep:ygrid[end]

  if length(yrange) < size(ygrid, 1)
    yrange = ygrid[1, t]:ystep:(ygrid[end, t] + 0.0001)
  end

  CoordInterpGrid((wrange, yrange), V[:, :, t], BCnearest, InterpLinear)

end

#######################################################################################
