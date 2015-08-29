################################################################################
###############################    SIMULATION    ###############################
################################################################################

function sim{T<:AbstractFloat}(wp::Array{T,6}, xgrid::Array{T,2},
  hgrid::Array{T,2}, agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1},
  wgrid_R::Array{T,2}, hgrid_R::Array{T,2}, ygrid_R::Array{T,1},
  yit::Array{T,2}, s_f_i::Array{T,3}, wp_R::Array{T,4}, pension::Array{T,1},
  r::T, λ::T)

  @printf "7. Simulate Consumption and Wealth Distribution\n"
  tW = size(yit,2); tR = size(wp_R,4)
  c_t = Array(Float64, (size(yit,1), tW+tR))
  w_t = similar(c_t); wp_t = similar(c_t); h_t = similar(c_t)

  # Initial conditions
  w_0 = 0.0; h_0 = mean(hgrid[:, 1]) + w_0/10 + rand(size(yit,1))
  w_t[:, 1] = w_0; h_t[:, 1] = h_0

  negcons = 0
  for t = 1:tW
    # INTERPOLATION
    wp_int =
      interpolateV(wp[:,:,:,:,:,t], xgrid[:,t], hgrid[:,t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]; ht = h_t[i, t]; yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, ht, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]
      h_t[i, t+1] = min((1-λ)*ht + λ*c_t[i, t], 0.005)
      c_t[i, t] < 0.0 ? negcons += 1 : 0
    end
  end
  negcons == 0 || println("\t$negcons negative consumption choices!")

  negcons = 0
  for t = (tW+1):(tW+tR)
    wp_int =
      interpolateV(wp_R[:,:,:,t-tW], wgrid_R[:,t-tW], hgrid_R[:,t-tW], ygrid_R)

    for  i = 1:size(yit,1)
      wt = w_t[i, t]; ht = h_t[i, t]
      xt = wt + pension[i]

      wp_t[i, t] = getValue(wp_int, [xt, ht, pension[i]])[1]
      c_t[i, t] = xt - wp_t[i, t]
      c_t[i, t] < 0.0 ? negcons += 1 : 0

      if t < tW + tR
        w_t[i, t+1] = r*wp_t[i, t]
        h_t[i, t+1] = min((1-λ)*ht + λ*c_t[i, t], 0.005)
      end
    end
  end
  negcons == 0 || println("\t$negcons negative c choices in retirement!")

  return c_t, h_t, w_t, wp_t
end
