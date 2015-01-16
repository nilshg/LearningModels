################################################################################
################################### SIMULATION #################################
################################################################################

function simulate(wp::Array, wgrid::Array, agrid::Array, bgrid::Array,
                  zgrid::Array, yit::Array, ymedian::Float64, s_f_i::Array,
                  wp_R::Array, wgrid_R::Array, ygrid_R::Array, pension::Array,
                  r::Float64)

  @printf "7. Simulate Consumption and Wealth Distribution\n"
  w_0 = 0.0
  c_t = Array(Float64, (agents*bs, tW+tR))
  w_t = similar(c_t)
  wp_t = similar(c_t)
  w_t[:, 1] = w_0
  tW = size(yit,2)
  tR = size(wgrid_R,2)

  @printf "\tSimulating %d periods of working life...\n" tW
  for t = 1:tW
    # INTERPOLATION
    wp_int = interpolatev(wp, wgrid, agrid, bgrid, zgrid, t)

    negconscounter = 0
    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]
      yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt
      wp_t[i, t] = wp_int[xt, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      if c_t[i, t] < 0.0 && negconscounter < 1
        @printf "\tWARNING: c=%.2f at wt=%.3f, yt=%.3f, t=%d\n" c_t[i,t] wt yt t
        negconscounter +=  1
      end
    end
  end

  @printf "\tSimulating %d periods of retirement...\n" tR
  for t = (tW+1):(tW+tR)
    negconscounter = 0
    wp_int = interpolatev(wp_R, wgrid_R, ygrid_R, t-tW)

    for  i = 1:size(yit,1)
      wt = w_t[i, t]
      yt = pension[i]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, yt]
      c_t[i, t] = xt - wp_t[i, t]

      if c_t[i, t] < 0.0 && negconscounter < 1
        @printf "\tWARNING: c=%.2f at wt=%.3f, yt=%.3f, t=%d\n" c_t[i,t] wt yt t
        negconscounter +=  1
      end

      if  t < tW + tR
        w_t[i, t+1] = r*wp_t[i, t]
      end
    end
  end

  return c_t, w_t, wp_t
end
