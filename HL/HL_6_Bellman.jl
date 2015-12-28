################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife{T<:AbstractFloat}(v::Array{T,6}, wp::Array{T,6},
  xgrid::Array{T,2}, hgrid::Array{T,2}, agrid::Array{T,1}, bgrid::Array{T,1},
  zgrid::Array{T,1}, stdy::Array{T,1}, k::Array{T,2}, r::T, δ::T, λ::T, ρ::T,
  c_over_x::Array{T,6})

  println("Solving for decision rules on $(prod(size(v)[1:4])) points on $(nprocs()) cores")

  wpnow = SharedArray(Float64, size(v)[1:5], pids=procs())
  vnow = SharedArray(Float64, size(v)[1:5], pids=procs())
  cxnow = SharedArray(Float64, size(v)[1:5], pids=procs())

  for t = (size(xgrid,2)-1):-1:1
    # INTERPOLATION
    v_interpol = interpolateV(v[:,:,:,:,:,t+1], xgrid[:,t+1], hgrid[:,t+1],
                                agrid, bgrid, zgrid)

    # MAXIMIZATION
    wmin = xgrid[1, t+1]
    @inbounds @sync @parallel for x = 1:size(xgrid,1)
      xt = xgrid[x, t]
      for h = 1:size(hgrid,1), a = 1:size(agrid,1), b = 1:size(bgrid,1)
        for z = 1:size(zgrid,1)
          ht = hgrid[h, t]; at = agrid[a]; bt = bgrid[b]; zt = zgrid[z]
          yn = Normal(at + bt*(t+1) + ρ*zt, stdy[t])

          (wpnow[x,h,a,b,z], vnow[x,h,a,b,z]) =
            bellOpt(xt, ht, at, bt, zt, wmin, v_interpol, yn, k[:,t], ρ, r, δ,λ)

          cxnow[x,h,a,b,z] = (xt-wpnow[x,h,a,b,z])/xt
        end
      end
    end
    v[:,:,:,:,:,t] = sdata(vnow); wp[:,:,:,:,:,t] = sdata(wpnow);
    c_over_x[:,:,:,:,:,t] = sdata(cxnow)
    mv = checkmonotonicity(v[:,:,:,:,:,t])
    println("\tMonotonicity violations: [w,h,a,b,z]=$(mv), t=$t")
  end

  return v, wp, c_over_x
end
