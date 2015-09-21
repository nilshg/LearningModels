################################################################################
################################### SIMULATION #################################
################################################################################

function sim{T<:AbstractFloat}(wp::Array{T,5}, xgrid::Array{T,2},
  agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1}, wgrid_R::Array{T,1},
  yit::Array{T,2}, s_f_i::Array{T,3}, pension::Array{T,1}, r::T, δ::T, σ::T,
  tR::Int64)

  @printf "Simulate Consumption and Wealth Distribution\n"
  tW = size(yit,2);
  c_t = Array(Float64, (size(yit,1), tW+tR))
  w_t = similar(c_t); wp_t = similar(c_t)

  # Initial wealth
  w_0 = 0.; w_t[:, 1] = w_0

  negcons = 0
  for t = 1:tW
    # INTERPOLATION
    wp_int = interpolateV(wp[:,:,:,:,t], xgrid[:,t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]; yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt

      wp_t[i, t] = getValue(wp_int, [xt, at, bt, zt])[1]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      c_t[i, t] < 0.0 ? negcons +=  1 : 0
    end
  end
  negcons == 0 || println("\t$negcons negative consumption choices!")

  function get_c_1{T<:AbstractFloat}(r::T, δ::T, x::T, y::T, σ::T, tR::Int64)
    numerator = 1 - 1/r*(r*δ)^(1/σ)
    denominator = 1 - (1/r*(r*δ)^(1/σ))^tR
    margprop = numerator/denominator
    pdvresources = y*((1-1/r^tR)/(1-1/r)) + x
    return margprop * pdvresources
  end

  function simulate_i{T<:AbstractFloat}(x::T, y::T, δ=δ, σ=σ, r=r, tR=tR)
    c_R_i = Array(Float64, tR); x_t = similar(c_R_i)
    c_R_i[1] = get_c_1(r, δ, x, y, σ, tR)[1]
    x_t[1] = x + y

    for t = 1:tR-1
      c_R_i[t+1] = (r*δ)^(1/σ)*c_R_i[t]
      x_t[t+1] = (x_t[t] - c_R_i[t])*r + y
    end
    return c_R_i, x_t
  end

  for i = 1:size(yit,1)
    (c_t[i, 41:70], w_t[i, 41:70]) =
      simulate_i(w_t[i, tW+1], pension[i])
  end

  return c_t, w_t, wp_t
end
