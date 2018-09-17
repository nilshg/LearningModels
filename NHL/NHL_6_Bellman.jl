################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array{T,4}, wp::Array{T,5}, xgrid::Array{T,2},
  agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1},
  stdy::Array{T,1}, k::Array{T,2}, r::T, δ::T, ρ::T, g_t::Array{T,1},
  σ::T, ξ::T) where T<:AbstractFloat

  wpnow = SharedArray{Float64}(size(v)[1:4], pids=procs());
  vnow = SharedArray{Float64}(size(v)[1:4], pids=procs());

  for t = (size(xgrid,2)-1):-1:1
    # INTERPOLATION
    v_interpol = interpolateV(v, xgrid[:,t+1], agrid, bgrid, zgrid)

    # MAXIMIZATION
    wmin = xgrid[1, t+1]
    @inbounds @sync @distributed for x = 1:size(xgrid,1)
      for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
        at = agrid[a]; bt = bgrid[b]; zt = zgrid[z]; xt = xgrid[x, t]
        yn = Normal(g_t[t+1] + at + bt*(t+1) + ρ*zt, stdy[t])

        (wpnow[x, a, b, z], vnow[x, a, b, z]) =
          bellOpt(xt, at, bt, zt, wmin, v_interpol, yn, k[:, t], ρ, r, δ, σ, ξ)

      end
    end
    wp[:,:,:,:,t] = sdata(wpnow); v = sdata(vnow)
  end
  return v, wp
end
