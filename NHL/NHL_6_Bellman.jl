################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array{Float64, 5}, wp::Array{Float64, 5},
                          c_over_x::Array{Float64, 5},
                          wgrid::Array{Float64, 2}, agrid::Array{Float64, 1},
                          bgrid::Array{Float64, 1}, zgrid::Array{Float64, 1},
                          stdy::Array, r::Float64, δ::Float64,
                          tW::Int64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving the problem on %d points\n" wpoints*apoints*bpoints*zpoints

  for t = (tW-1):-1:1
    @printf "    Period %d/%d, w = [%.2f, %.2f]\n" t tW wgrid[1,t] wgrid[end,t]

    # INTERPOLATION
    v_interpol = interpolatev(v, wgrid, agrid, bgrid, zgrid, t+1)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    for w = 1:size(wgrid, 1)
      for a = 1:size(agrid, 1)
        for b = 1:size(bgrid, 1)
          for z = 1:size(bgrid, 1)
            size(agrid,2)==1 ? at = agrid[a] : at = agrid[a,t+1]
            size(bgrid,2)==1 ? bt = bgrid[b] : bt = bgrid[b,t+1]
            size(zgrid,2)==1 ? zt = zgrid[z] : zt = zgrid[z,t+1]
            wt = wgrid[w, t]
            yt = exp(at + bt*t + zt)
            yln = LogNormal(at + bt*(t+1) + zt, stdy[t])

            (wt + yt - wmin/r > 0.01) || error("Cash on hand is wt + yt, BC is wmin/r")

            (wpopt, vopt) = bellOpt(wt, yt, at, bt, zt, wmin,
                                    v_interpol, yln, r, δ)

            wpopt < wt + yt || @printf "\tw'>xt, w=%d,a=%d,b=%d,z=%d,t=%d\n" w a b z t

            c = wgrid[w, t] + yt - wpopt
            c_over_x[w, a, b, z, t] = c/(wgrid[w, t] + yt - wmin/r)

            wp[w, a, b, z, t] = wpopt
            v[w, a, b, z, t] = vopt
          end
        end
      end
    end
  end

  return v, wp, c_over_x
end
