################################################################################
###############################    SIMULATION    ###############################
################################################################################

function simulate(wp::Array{Float64, 6}, wgrid::Array{Float64},
                  hgrid::Array{Float64}, agrid::Array{Float64},
                  bgrid::Array{Float64}, zgrid::Array{Float64},
                  wgrid_R::Array{Float64}, hgrid_R::Array{Float64},
                  ygrid_R::Array{Float64}, yit::Array{Float64, 2},
                  ymedian::Float64, s_f_i::Array{Float64, 3},
                  wp_R::Array, r::Float64, λ::Float64)

  @printf "7. Simulate Consumption and Wealth Distribution\n"

  tW = size(yit,2)
  tR = size(wp_R,4)

  w_0 = 0
  h_0 = mean(hgrid[:, 1]) + w_0/10 + rand(size(yit,1))

  c_t = Array(Float64, (size(yit,1), tW+tR))
  w_t = similar(c_t)
  wp_t = similar(c_t)
  h_t = similar(c_t)

  w_t[:, 1] = w_0
  h_t[:, 1] = h_0

  @printf "\tSimulating %d periods of working life...\n" tW
  for t = 1:tW
    # INTERPOLATION
    wp_int = interpolatev(wp, wgrid, hgrid, agrid, bgrid, zgrid, t)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]
      ht = h_t[i, t]
      yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, ht, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]

      if t < tW
        w_t[i, t+1] = r*wp_t[i, t]
        h_t[i, t+1] = (1-λ)*ht + λ*c_t[i, t]
      end
    end
  end

  pension = zeros(size(yit,1), 1)

  for i = 1:size(yit,1)
    yt = yit[i, tW]

    if yt < 0.3*ymedian
      yfixed = 0.9*yt
    elseif yt <= 2.0*ymedian
      yfixed = 0.27*ymedian + 0.32(yt - 0.3*ymedian)
    elseif yt <= 4.1*ymedian
      yfixed = 0.81*ymedian + 0.15*(yt - 2*ymedian)
    elseif yt > 4.1*ymedian
      yfixed = 1.1*yt
    end
    pension[i] = 0.715*yfixed
  end

  @printf "\t\Simulating %d periods of retirement...\n" tR

  for t = (tW+1):(tW+tR)
    negconscounter = 0
    wp_int = interpolatev(wp_R, wgrid_R, hgrid_R, ygrid_R, t-tW)

    for  i = 1:size(yit,1)
      wt = w_t[i, t]
      ht = h_t[i, t]
      yt = pension[i][1]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, ht, yt]
      c_t[i, t] = xt - wp_t[i, t]

      if c_t[i, t] < 0.0 && negconscounter < 5
        @printf "Warning, neg cons at wt=%.3f, ht=%.3f, yt=%.3f, t=%d\n" wt ht yt t
        negconscounter +=  1
      end

      if t < tW + tR
        w_t[i, t+1] = r*wp_t[i, t]
        h_t[i, t+1] = (1-λ)*ht + λ*c_t[i, t]
      end
    end
  end

  return c_t, h_t, w_t, wp_t
end
