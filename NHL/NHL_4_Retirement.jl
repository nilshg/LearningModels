###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################

function solveRetirement(wgrid_R::Array, ygrid_R::Array, r::Float64, δ::Float64)

  @printf "4. Solving the retirement problem\n"
  tic()
  tR = size(wgrid_R,2)
  v_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))
  wp_R = similar(v_R)

  # Value of last period of retirement
  for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid,1)
      wp_R[w, y, tR] = 0.0
      v_R[w, y, tR] = u(wgrid_R[w, tR] + ygrid_R[y])
    end
  end

  # Compute the period tR-1 solution exactly:
  wmin = wgrid_R[1, tR]
  for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid,1)
      wt = wgrid_R[w, tR-1]
      yt = ygrid_R[y]
      xt = wt + yt

      Blmn(wp::Float64) = -(u(xt - wp) + δ*u(r*wp+yt))

      optimum = optimize(Blmn, wmin/r, xt)

      v_R[w, y, tR-1] = -(optimum.f_minimum)
      wp_R[w, y, tR-1] = optimum.minimum
    end
  end

  for t = (tR-2):-1:1
    v_R_interpol = interpolatev(v_R, wgrid_R, ygrid_R, t+1)
    wmin = wgrid_R[1, t+1]
    for w = 1:size(wgrid_R,1)
     for y = 1:size(ygrid,1)
        (wp_R[w, y, t], v_R[w, y, t]) =
          bellOpt_R(wgrid_R[w, t], ygrid_R[y], wmin, v_R_interpol, r, δ)
      end
    end
  end
  @printf "\tSolving retirement problem took %.1f seconds.\n" toq()
  return v_R, wp_R
end
