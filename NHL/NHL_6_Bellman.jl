################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife{T<:AbstractFloat}(v::Array{T,5}, wp::Array{T,5},
  xgrid::Array{T,2}, agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1},
  stdy::Array{T,1}, k::Array{T,2}, r::T, δ::T, ρ::T, c_over_x::Array{T,5})

  println("Solve for policy function on $(prod(size(v)[1:4])) points")
  wpnow = SharedArray(Float64, size(v)[1:4], pids=procs()); vnow=similar(wpnow)
  cxnow = similar(wpnow)

  for t = (size(xgrid,2)-1):-1:1
    # INTERPOLATION
    v_interpol = interpolateV(v[:,:,:,:,t+1], xgrid[:,t+1], agrid, bgrid, zgrid)

    # MAXIMIZATION
    wmin = xgrid[1, t+1]
    @inbounds @sync @parallel for x = 1:size(xgrid,1)
      for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
        at = agrid[a]; bt = bgrid[b]; zt = zgrid[z]; xt = xgrid[x, t]
        yn = Normal(at + bt*(t+1) + ρ*zt, stdy[t])

        (wpnow[x, a, b, z], vnow[x, a, b, z]) =
          bellOpt(xt, at, bt, zt, wmin, v_interpol, yn, k[:, t], ρ, r, δ)

        cxnow[x,a,b,z] = (xt-wpnow[x,a,b,z])/xt
      end
    end
    wp[:,:,:,:,t] = sdata(wpnow); v[:,:,:,:,t] = sdata(vnow)
    c_over_x[:,:,:,:,t] = sdata(cxnow)
  end
  return v, wp, c_over_x
end
