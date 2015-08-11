################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array, wp::Array, wgrid::Array, hgrid::Array,
                          agrid::Array, bgrid::Array, zgrid::Array,
                          stdy::Array, r::Float64, k::Array, δ::Float64,
                          λ::Float64, ρ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving the problem on %d points\n" length(v)/size(wgrid,2)

  for t = (size(wgrid,2)-1):-1:1
    @printf "     Period: %d out of %d\n" t size(wgrid,2)

    # INTERPOLATION
    v_interpol = interpolateV(v[:,:,:,:,:,t+1],
                              wgrid[:,t+1], hgrid[:,t+1], agrid, bgrid, zgrid)

    wpnow = SharedArray(Float64, (size(v[:,:,:,:,:,1])), pids=procs())
    vnow = similar(wpnow)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    tic()
    @inbounds @sync @parallel for w = 1:size(wgrid,1)
      wt = wgrid[w, t]
      for h = 1:size(hgrid,1), a = 1:size(agrid,1), b = 1:size(bgrid,1)
        for z = 1:size(zgrid,1)
          ht = hgrid[h, t]
          at = agrid[a]
          bt = bgrid[b]
          zt = zgrid[z]
          yt = exp(at + bt*t + zt)
          yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

          if wt + yt - wmin/r < 0.01
            error("Error: Cash on hand is too low, at w, a, b, z")
          end

          (wpopt, vopt) = bellOpt(wt, ht, yt, at, bt, zt, wmin,
                                  v_interpol, yln, k[:,t], ρ, r, δ, λ)

          wpnow[w,h,a,b,z] = wpopt
          vnow[w,h,a,b,z] = vopt

          wpopt < wt + yt || @printf "NC @ w=%d,h=%d,y=%d,t=%d" w h yt t
        end
      end
    end
    v[:,:,:,:,:,t] = sdata(vnow)
    wp[:,:,:,:,:,t] = sdata(wpnow)

    mv = checkmonotonicity(v[:,:,:,:,:,t])
    @printf "\tMonotonicity violations: w=%d, h=%d, a=%d, b=%d, z=%d\n" mv[1] mv[2] mv[3] mv[4] mv[5]
  end
  @printf "\tMaximization took %d seconds\n" toq()

  return v, wp
end
