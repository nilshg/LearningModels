###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################
@printf "5. Solving the problem for the last period of work\n"

wp = zeros(wgridpoints, hgridpoints, agridpoints, bgridpoints, zgridpoints, T)
V = zeros(wgridpoints, hgridpoints, agridpoints, bgridpoints, zgridpoints, T)

valueRETIRE = interpolateV_R(V_R, wgrid_R, hgrid_R, ygrid_R, 1)

for b = 1:bgridpoints
  tic()
  for a = 1:agridpoints
    for z = 1:zgridpoints

      at = agrid[1, T]
      bt = bgrid[1, T]
      zt = zgrid[1, T]

      ytemp = exp(at + bt*T + zt) + 0.4

      if ytemp < 0.3*ymedian
        yfixed = 0.9*ytemp
      elseif ytemp <= 2.0*ymedian
        yfixed = 0.27*ymedian + 0.32(ytemp - 0.3*ymedian)
      elseif ytemp < 4.1*ymedian
        yfixed = 1.1*ymedian
      end
      yfixed = 0.715*yfixed

      wmin = wgrid_R[1, 1]

      for w = 1: wgridpoints
        for h = 1:hgridpoints
          for y = 1:ygridpoints

          wt = wgrid[w, T]
          ht = hgrid[h, T]
          yt = ytemp - 0.4

          (wpopt, Vopt) = BellOpt_TRANS(wt, ht, yt, yfixed, wgrid_R[1, 1],
                                        valueRETIRE, u, R, dbeta, lambda, 40)

          V[w, h, a, b, z, T] = Vopt
          wp[w, h, a, b, z, T] = wpopt
          end
        end
      end
    end
  end
  @printf "\tβ = %d/%d, one β takes %.2f seconds\n" b bgridpoints toq()
end
