###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################

function solveTransition(v_R::Array, wgrid_R::Array, hgrid_R::Array, ygrid_R::Array,
                         wgrid::Array, hgrid::Array, agrid::Array, bgrid::Array,
                         zgrid::Array, ymedian::Float64, r::Float64, δ::Float64,
                         λ::Float64)

  @printf "5. Solving the problem for the last period of work\n"
  tic()
  tW = size(wgrid,2)
  wp = Array(Float64,
             size(wgrid,1), size(hgrid,1), size(agrid,1), size(bgrid,1),
             size(zgrid,1), tW)
  v = similar(wp)

  # INTERPOLATION
  valueRETIRE = interpolateV(v_R[:,:,:,1], wgrid_R[:,1], hgrid_R[:,1], ygrid_R)

  # MAXIMIZATION
  wmin = wgrid_R[1, 1]

  for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
    at = agrid[a]
    bt = bgrid[a]
    zt = 0

    yt = exp(at + bt*tW + zt)

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

    for w = 1: size(wgrid,1), h = 1:size(hgrid,1)
      wt = wgrid[w, tW]
      ht = hgrid[h, tW]

      if wt + yfixed < wmin
        @printf "\tWarning: Cash in hand too low at w=%d, a=%d, b=%d, z=%d\n" w a b z
      end

      (wpopt, vopt) = bellOpt_TRANS(wt, ht, yt, yfixed, wmin, valueRETIRE, r, δ, λ)

      v[w, h, a, b, z, tW] = vopt
      wp[w, h, a, b, z, tW] = wpopt

      if isnan(vopt)
        @printf "Warning: V is NaN at w=%d, h=%d, y=%d" wt ht yt
      end
      if wpopt > wt + yt
        @printf "Warning: Negative c at w=%d, h=%d, y=%d" wt ht yt
      end
    end
  end
  mv = checkmonotonicity(v[:,:,:,:,:,end])
  @printf "\tMonotonicity violations: w=%d, h=%d, a=%d, b=%d, z=%d\n" mv[1] mv[2] mv[3] mv[4] mv[5]
  @printf "\tTransition period problem solved in %.1f seconds.\n" toq()

  return v, wp
end
