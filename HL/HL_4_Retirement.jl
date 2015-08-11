###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################

function solveRetirement(wgrid_R::Array, hgrid_R::Array, ygrid_R::Array,
                         r::Float64, δ::Float64, λ::Float64)

  @printf "4. Solving the retirement problem\n"
  tic()
  tR = size(wgrid_R,2)
  v_R = Array(Float64, (size(wgrid_R,1), size(hgrid_R,1), size(ygrid_R,1), tR))
  wp_R = similar(v_R)

  # Value of last period of retirement
  for w = 1:size(wgrid_R,1), h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
    wt = wgrid_R[w, tR]
    ht = hgrid_R[h, tR]
    yt = ygrid_R[y]
    xt = wt + yt
    v_R[w, h, y, tR] = u(xt, ht)
    wp_R[w, h, y, tR] = 0.0
  end

  # Compute the period tR-1 solution exactly:
  wmin = wgrid_R[1, tR]
  for w = 1:size(wgrid_R,1), h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
    wt = wgrid_R[w, tR-1]
    ht = hgrid_R[h, tR-1]
    yt = ygrid_R[y]
    xt = wt + yt

    Blmn(wp::Float64) = -(u(xt - wp/r, ht)
                          + δ*u(wp+yt, (1-λ)*h + λ*(xt-wp/r)))

    optimum = optimize(Blmn, wmin, xt-0.01)

    v_R[w, h, y, tR-1] = -(optimum.f_minimum)
    wp_R[w, h, y, tR-1] = optimum.minimum
  end

  mv = checkmonotonicity(v_R[:,:,:,tR-1])
  @printf "\tMonotonicity violations: w=%d, h=%d, y=%d, t=%d\n" mv[1] mv[2] mv[3] tR-1

  for t = (tR-2):-1:1
    # INTERPOLATION
    v_R_interpol = interpolateV(v_R[:,:,:,t+1], wgrid_R[:,t], hgrid_R[:,t], ygrid_R)

    # MAXIMIZATION
    neg_c_counter = zeros(30)
    nan_counter = zeros(30)
    wmin = wgrid_R[1, t+1]
    for w = 1:size(wgrid_R,1), h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
      wt = wgrid_R[w, t]
      ht = hgrid_R[h, t]
      yt = ygrid_R[y]

      wt + yt > wmin/r ? 0 : error("xt < wmin/r!")

      (wpopt, vopt) = bellOpt_R(wt, ht, yt, wmin, v_R_interpol, r, δ, λ)

      v_R[w, h, y, t] = vopt
      wp_R[w, h, y, t] = wpopt

      if isnan(vopt)
        nan_counter[t] += 1
      end
      if wpopt > wt + yt
        neg_c_counter[t] += 1
      end
    end

    mv = checkmonotonicity(v_R[:,:,:,t])
    @printf "\tMonotonicity violations: w=%d, h=%d, y=%d, t=%d\n" mv[1] mv[2] mv[3] t
  end
  @printf "\tSolving retirement problem took %.1f seconds.\n" toq()

  return v_R, wp_R
end
