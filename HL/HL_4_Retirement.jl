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
  for w = 1:size(wgrid_R,1)
    for h = 1:size(hgrid_R,1)
      for y = 1:size(ygrid_R,1)
        wt = wgrid_R[w, tR]
        ht = hgrid_R[h, tR]
        yt = ygrid_R[y]
        xt = wt + yt
        v_R[w, h, y, tR] = u(xt, ht)
        wp_R[w, h, y, tR] = 0.0
      end
    end
  end

  # Compute the period tR-1 solution exactly:
  wmin = wgrid_R[1, tR]
  for w = 1:size(wgrid_R,1)
    for h = 1:size(hgrid_R,1)
      for y = 1:size(ygrid_R,1)
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
    end
  end

  (dim1, dim2, dim3) = checkmonotonicity(v_R, tR-1)
  @printf "\tMonotonicity violations: w=%d, h=%d, y=%d, t=%d\n" dim1 dim2 dim3 tR-1

  for t = (tR-2):-1:1
    # INTERPOLATION
    v_R_interpol = interpolatev(v_R, wgrid_R, hgrid_R, ygrid_R, t+1)

    # MAXIMIZATION
    wmin = wgrid_R[1, t+1]
    for w = 1:size(wgrid_R,1)
      for h = 1:size(hgrid_R,1)
        for y = 1:size(ygrid_R,1)
          wt = wgrid_R[w, t]
          ht = hgrid_R[h, t]
          yt = ygrid_R[y]

          (wpopt, vopt) = bellOpt_R(wt, ht, yt, wmin, v_R_interpol, u, r, δ, λ)

          v_R[w, h, y, t] = vopt
          wp_R[w, h, y, t] = wpopt

          if isnan(vopt)
            @printf "Warning: V is NaN at w=%d, h=%d, y=%d, t=%d" w h y t
          end
          if wpopt > wt + yt
            @printf "Warning: Negative c at w=%d, h=%d, y=%d, t=%d" w h y t
          end
        end
      end
    end
    (dim1, dim2, dim3) = checkmonotonicity(v_R, t)
    @printf "\tMonotonicity violations: w=%d, h=%d, y=%d, t=%d\n" dim1 dim2 dim3 t
  end
  @printf "\tSolving retirement problem took %.1f seconds.\n" toq()

  return v_R, wp_R
end
