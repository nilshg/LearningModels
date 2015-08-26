################################################################################
################################### SIMULATION #################################
################################################################################

function sim(wp::Array, wgrid::Array, agrid::Array, bgrid::Array, zgrid::Array,
  yit::Array, ymedian::Float64, s_f_i::Array, wp_R::Array, wgrid_R::Array,
  ygrid_R::Array, pension::Array, r::Float64)

  @printf "7. Simulate Consumption and Wealth Distribution\n"
  w_0 = 0.0; tW = size(yit,2); tR = size(wgrid_R,2)
  c_t = Array(Float64, (size(yit,1), tW+tR))
  w_t = similar(c_t); wp_t = similar(c_t)
  w_t[:, 1] = w_0

  negconscounter = 0
  for t = 1:tW
    # INTERPOLATION
    wp_int =
      interpolateV(wp[:, :, :, :, t], wgrid[:, t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]
      yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt
      wp_t[i, t] = getValue(wp_int, [xt, at, bt, zt])[1]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      if c_t[i, t] < 0.0
        negconscounter +=  1
      end
    end
  end
  negconscounter == 0 || @printf "\t%d neg. c choices!\n" negconscounter

  function get_c_1{T<:AbstractFloat}(r::T, δ::T, x::T, y::T, σ::T, tR::Int64)
    numerator = 1 - 1/r*(r*δ)^(1/σ)
    denominator = 1 - (1/r*(r*δ)^(1/σ))^tR
    margprop = numerator/denominator
    pdvresources = y*((1-1/r^tR)/(1-1/r)) + x
    return margprop * pdvresources
  end

  function simulate_i{T<:AbstractFloat}(x::T, y::T, δ::T, σ::T, r::T, tR::Int64)
    c_R_i = Array(Float64, tR)
    x_t = similar(c_R_i)
    c_R_i[1] = get_c_1(r, δ, x, y, σ, tR)[1]
    x_t[1] = x + y

    for t = 1:tR-1
      c_R_i[t+1] = (r*δ)^(1/σ)*c_R_i[t]
      x_t[t+1] = (x_t[t] - c_R_i[t])*r + y
    end
    return c_R_i, x_t
  end

  for i = 1:size(yit,1)
    (a,b) =
      simulate_i(w_t[i, tW+1], pension[i], δ, σ, r, tR)
    c_t[i, 41:70] = a
    w_t[i, 41:70] = b
  end

  return c_t, w_t, wp_t
end
