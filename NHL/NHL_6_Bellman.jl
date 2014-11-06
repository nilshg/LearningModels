###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################
@printf "6. Solving the problem for the last period of work\n"
tic()
wp = zeros(wgridpoints, agridpoints, bgridpoints, zgridpoints, T)
V = zeros(wgridpoints, agridpoints, bgridpoints, zgridpoints, T)

# INTERPOLATION
valueRETIRE = interpolateV(V_R, wgrid_R, ygrid_R, 1)

# MAXIMIZATION
wmin = wgrid_R[1, 1]

for a = 1:agridpoints
  for b = 1:bgridpoints
    for z = 1:zgridpoints
      at = agrid[a, T]
      bt = bgrid[b, T]
      zt = 0

      yt = exp(at + T*bt + zt)

      if yt < 0.3*ymedian
        yfixed = 0.9*yt
      elseif yt <= 2.0*ymedian
        yfixed = 0.27*ymedian + 0.32(yt - 0.3*ymedian)
      elseif yt <= 4.1*ymedian
        yfixed = 0.81*ymedian + 0.15*(yt - 2*ymedian)
      elseif yt > 4.1*ymedian
        yfixed = 1.1*yt
      end
      yfixed = 0.715*yfixed

      for w = 1: wgridpoints
        wt = wgrid[w, T]

        (wpopt, Vopt) = BellOpt_TRANS(wt, yt, yfixed, wmin, valueRETIRE, u, R, dbeta, T)

        V[w, a, b, z, T] = Vopt
        wp[w, a, b, z, T] = wpopt
      end
    end
  end
end
@printf "\tTransition period problem solved in %.1f seconds.\n" toq()
