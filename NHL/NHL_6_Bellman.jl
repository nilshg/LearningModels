################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array{Float64, 5}, wp::Array{Float64, 5},
           c_over_x::Array{Float64, 5}, wgrid::Array{Float64, 2},
           agrid::Array{Float64}, bgrid::Array{Float64}, zgrid::Array{Float64},
           stdy::Array{Float64}, k::Array, r::Float64, δ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving on %d points..." length(v[:, :, :, :, end])

  for t = (size(wgrid,2)-1):-1:1
    mod(t,10) != 0 || @printf "%d/%d, " t size(wgrid,2)

    # INTERPOLATION
    v_interpol =
      interpolateV(v[:, :, :, :, t+1], wgrid[:, t+1], agrid, bgrid, zgrid)

    wpnow = SharedArray(Float64, (size(v[:, :, :, :, 1])), pids=procs())
    cxnow = similar(wpnow); vnow = similar(wpnow)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    @inbounds @sync @parallel for w = 1:size(wgrid,1)
      for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
        size(agrid,2)==1 ? at = agrid[a] : at = agrid[a, t]
        size(bgrid,2)==1 ? bt = bgrid[b] : bt = bgrid[b, t]
        size(zgrid,2)==1 ? zt = zgrid[z] : zt = zgrid[z, t]
        wt = wgrid[w, t]
        yt = exp(at + bt*t + zt)
        yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

        (wt + yt - wmin/r > 0.01) || error("Cash on hand is $wt + $yt, BC is $wmin/r")

        (wpopt, vopt) =
          bellOpt(wt, yt, at, bt, zt, wmin, v_interpol, yln, k[:, t], ρ, r, δ)

        wpopt < wt + yt || @printf "\tw'>xt, w=%d,a=%d,b=%d,z=%d,t=%d\n" w a b z t

        c = wgrid[w, t] + yt - wpopt
        cxnow[w, a, b, z] = c/(wgrid[w, t] + yt - wmin/r)

        wpnow[w, a, b, z] = wpopt
        vnow[w, a, b, z] = vopt
      end
    end
    wp[:, :, :, :, t] = sdata(wpnow)
    v[:, :, :, :, t] = sdata(vnow)
    c_over_x[:, :, :, :, t] = sdata(cxnow)
  end

  return v, wp, c_over_x
end
