################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################

function solveWorkingLife(v::Array, wp::Array, wgrid::Array, hgrid::Array,
                          agrid::Array, bgrid::Array, zgrid::Array, stdy::Array,
                          r::Float64, δ::Float64, λ::Float64, ρ::Float64)

  @printf "6. Recursively solve for optimal decision rules\n"
  @printf "\tSolving the problem on %d points\n" length(v)/size(wgrid,2)

  for t = (size(wgrid,2)-1):-1:1
    @printf "     Period: %d out of %d\n" t size(wgrid,2)
    @printf "\tw = [%.2f, %.2f], h = [%.2f, %.2f]\n" wgrid[1, t] wgrid[end, t] hgrid[1, t] hgrid[end, t]

    # INTERPOLATION
    v_interpol = interpolatev(v, wgrid, hgrid, agrid, bgrid, zgrid, t+1)

    # MAXIMIZATION
    wmin = wgrid[1, t+1]
    tic()
    isnancounter = 0
    for a = 1:size(agrid,1)
      for b = 1:size(bgrid,1)
        for z = 1:size(zgrid,1)
          at = agrid[a]
          bt = bgrid[b]
          zt = zgrid[z]
          yt = exp(at + bt*t + zt)
          yln = LogNormal(at + bt*(t+1) + ρ*zt, stdy[t])

          for h = 1:hpoints
            for w = 1:wpoints
              ht = hgrid[h, t]
              wt = wgrid[w, t]

              if wt + yt - wmin/r < 0.01
                error("Error: Cash on hand is too low, at w, a, b, z")
              end

              (wpopt, vopt) = bellOpt(wt, ht, yt, at, bt, zt, wmin,
                                      v_interpol, yln, r, δ, λ)

              wp[w, h, a, b, z, t] = wpopt
              v[w, h, a, b, z, t] = vopt

              wpopt < wt + yt || @printf "NC @ w=%d,h=%d,y=%d,t=%d" w h y t
            end
          end
        end
      end
    end
    (dim1, dim2, dim3, dim4, dim5) = checkmonotonicity(v, t)
    @printf "\tMonotonicity violations: w=%d, h=%d, a=%d, b=%d, z=%d\n" dim1 dim2 dim3 dim4 dim5
    @printf "\tMaximization took %d seconds\n" toq()
  end

  return v, wp
end
