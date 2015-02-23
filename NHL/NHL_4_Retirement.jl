###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################

function solveRetirement(wgrid_R::Array, ygrid_R::Array, r::Float64, δ::Float64)

  @printf "4. Solving the retirement problem\n"

  tR = size(wgrid_R,2)
  v_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))
  wp_R = similar(v_R)
  c_over_x = similar(v_R)

  # Value of last period of retirement
  for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid_R,1)
      wp_R[w, y, tR] = 0.0
      v_R[w, y, tR] = u(wgrid_R[w, tR] + ygrid_R[y])
      c_over_x[w, y, tR] = 1.0
    end
  end

  # Compute the period tR-1 solution exactly:
  wmin = wgrid_R[1, tR]
  for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid_R,1)
      wt = wgrid_R[w, tR-1]
      yt = ygrid_R[y]
      xt = wt + yt

      Blmn(wp::Float64) = -(u(xt - wp) + δ*u(r*wp+yt))

      optimum = optimize(Blmn, wmin/r, xt)

      v_R[w, y, tR-1] = -(optimum.f_minimum)
      wp_R[w, y, tR-1] = optimum.minimum

      c = xt - wp_R[w, y, tR-1]
      c_over_x[w, y, tR-1] = c/(xt - wmin/r)
    end
  end

  for t = (tR-2):-1:1
    v_R_interpol = interpolatev(v_R, wgrid_R, ygrid_R, t+1)
    wmin = wgrid_R[1, t+1]
    for w = 1:size(wgrid_R,1)
     for y = 1:size(ygrid_R,1)
        (wp_R[w, y, t], v_R[w, y, t]) =
          bellOpt_R(wgrid_R[w, t], ygrid_R[y], wmin, v_R_interpol, r, δ)

        xt = wgrid_R[w, t] + ygrid_R[y]
        c = xt - wp_R[w, y, t]
        c_over_x[w, y, t] = c/(xt - wmin/r)
      end
    end
  end

  return v_R, wp_R, c_over_x
end


################################################################################

function solveRetirement(wgrid_R::Array, ygrid_R::Array, r::Float64, δ::Float64,
                         σ::Float64)

  function get_c_1(r::Float64, δ::Float64, x::Float64, y::Float64,
                   σ::Float64, tR::Int64)
    Rinv = 1/r
    Rbsig = (r*δ)^(1/σ)
    numerator = 1 - Rinv*Rbsig
    denominator = 1 - (Rinv*Rbsig)^tR
    margprop = numerator/denominator
    pdvlabour = y*((1-Rinv^tR)/(1-Rinv))
    pdvresources = pdvlabour + x
    c_1 = margprop * pdvresources

    return c_1
  end

  function simulate(x::Float64, y::Float64, δ::Float64, σ::Float64,
                    r::Float64, tR::Int64)
    c_t = Array(Float64, tR)
    x_t = similar(c_t)
    sum_u = 0.0

    c_t[1] = get_c_1(r, δ, x, y, σ, tR)[1]
    x_t[1] = x + y

    for t = 1:tR-1
      c_t[t+1] = (r*δ)^(1/σ)*c_t[t]
      x_t[t+1] = (x_t[t] - c_t[t])*r + y
      sum_u += (δ^t)*(c_t[t]^(1-σ)/(1-σ))
    end
    sum_u += (δ^tR)*(c_t[tR]^(1-σ)/(1-σ))

    return x_t, sum_u
  end

  wp_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), size(wgrid_R,2)))
  v_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), 1))

  for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid_R,1)
      xt = wgrid_R[w]
      yt = ygrid_R[y]

      (wp_R[w, y, :], v_R[w, y]) =
        simulate(xt, yt, δ, σ, r, tR)
    end
  end
  return wp_R, v_R
end
