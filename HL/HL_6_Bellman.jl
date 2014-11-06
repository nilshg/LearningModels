################################################################################
############################## BELLMAN RECURSION ###############################
################################################################################
@printf "7. Recursively solve for optimal decision rules"

for t = (T-1):-1:1
  @printf "     Period: %d out of %d\n" t T
  @printf "\tw = [%.2f, %.2f], y = [%.2f, %.2f], h = [%.2f, %.2f]\n" wgrid[1, t] wgrid[end, t] minimum(Ybelief[:,1:Qt[t],t]) maximum(Ybelief[:,:,t]) hgrid[1, t] hgrid[end, t]

  # INTERPOLATION
  V_interpol = interpolateV(V, wgrid, hgrid, agrid, bgrid, zgrid, t+1)

  # MAXIMIZATION
  wmin = wgrid[1, t+1]

  for b = 1:bgridpoints
    tic()
    for a = 1:agridpoints
      for z = 1:zgridpoints
        for y = 1:ygridpoints
          for h = 1:hgridpoints
            for w = 1:wgridpoints
              ht = hgrid[h, t]
              wt = wgrid[w, t]
              at = agrid[a, t]
              bt = bgrid[b, t]
              zt = zgrid[z, t]
              yt = a + b*(t+1) + z

              Ydist = LogNormal(yt, stdy[t])

              if wt + yt - wmin/R < 0.01
                @printf "Error: Cash on hand is too low, at w=%d, h=%d, y=%d" w h y
              pause
              end

              (wpopt, Vopt) = BellOpt(wt, ht, yt, at, bt, zt, wmin,
                                      V_interpol, Ydist, u, R, lambda, dbeta, t)

              wp[w, h, a, b, z, t] = wpopt
              V[w, h, a, b, z, t] = Vopt
            end
          end
        end
      end
    end
    @printf "\tβ is %d/%d, one β takes %.2f seconds\n" b bgridpoints toq()
  end
end
