################################################################################
################################### SIMULATION #################################
################################################################################
using Distributions

function get_percentiles(w::Array{Float64,2})
  prct = zeros(90,4)
  for i = 1:90
    prct[i,1] = percentile(mean(w[:,6:15], dims = 2)[:], i)
    prct[i,2] = percentile(mean(w[:,16:25], dims = 2)[:], i)
    prct[i,3] = percentile(mean(w[:,26:35], dims = 2)[:], i)
    prct[i,4] = percentile(mean(w[:,6:35], dims = 2)[:], i)
  end
  return prct
end

function sim(wp::Array{T,5}, xgrid::Array{T,2}, agrid::Array{T,1},
  bgrid::Array{T,1}, zgrid::Array{T,1}, wgrid_R::Array{T,1},
  yit::Array{T,2}, s_f_i::Array{T,3}, pension::Array{T,1}, r::T, δ::T, σ::T,
  tR::Int64, p::AbstractString) where T<:AbstractFloat

  println("Simulate Consumption and Wealth Distribution")
  tW = size(yit,2);
  c_t = Array{Float64}(undef, size(yit,1), tW+tR)
  w_t = similar(c_t); wp_t = similar(c_t)

  # Initial wealth
  #w_0 = 0.; w_t[:, 1] = w_0
  include(p*"SCF_1983_age_distribution.jl")
  include(p*"SCF_percentiles.jl")

  for i = 1:size(w_t,1)
    w_t[i,1] = sample(w_0_weight[:,1], weights(w_0_weight[:,2]))
  end

  negcons = 0
  for t = 1:tW
    # INTERPOLATION
    wp_int = interpolateV(wp[:,:,:,:,t], xgrid[:,t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]; yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      c_t[i, t] < 0.0 ? negcons +=  1 : 0
    end
  end
  negcons == 0 || println("\t$negcons negative consumption choices!")

  function get_c_1(r, δ, x, y, σ, tR::Int64)
    numerator = 1 - 1/r*(r*δ)^(1/σ)
    denominator = 1 - (1/r*(r*δ)^(1/σ))^tR
    margprop = numerator/denominator
    pdvresources = y*((1-1/r^tR)/(1-1/r)) + x
    return margprop * pdvresources
  end

  function simulate_i(x, y, δ=δ, σ=σ, r=r, tR=tR)
    c_R_i = Array{Float64}(undef, tR); x_t = similar(c_R_i)
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

  avgy = mean(yit); w_t /= avgy; c_t /= avgy
  percentiles = get_percentiles(w_t)

  return c_t, w_t, wp_t, percentiles
end

function sim(wp::Array{T,5}, wp_R::Array{T,3}, xgrid::Array{T,2},
  agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1}, wgrid_R::Array{T,1},
  ygrid_R::Array{T,1}, yit::Array{T,2}, s_f_i::Array{T,3}, pension::Array{T,1},
  r::T, δ::T, σ::T, tR::Int64, p::AbstractString) where T<:AbstractFloat

  println("Simulate Consumption and Wealth Distribution")
  tW = size(yit,2);
  c_t = Array{Float64}(size(yit,1), tW+tR)
  w_t = similar(c_t); wp_t = similar(c_t)

  # Initial wealth
  #w_0 = 0.; w_t[:, 1] = w_0
  include(p*"SCF_1983_age_distribution.jl")
  include(p*"SCF_percentiles.jl")

  for i = 1:size(w_t,1)
    w_t[i,1] = sample(w_0_weight[:,1], weights(w_0_weight[:,2]))
  end

  negcons = 0
  for t = 1:tW
    # INTERPOLATION
    wp_int = interpolateV(wp[:,:,:,:,t], xgrid[:,t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wt = w_t[i, t]; yt = yit[i, t]
      (at, bt, zt) = s_f_i[:, i, t]
      xt = wt + yt

      wp_t[i, t] = wp_int[xt, at, bt, zt]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      c_t[i, t] < 0.0 ? negcons +=  1 : 0
    end
  end
  negcons == 0 || println("\t$negcons negative consumption choices!")

  for t = tW+1:tR
    # INTERPOLATION
    wp_int = interpolateV(wp[:,:,:,:,t], xgrid[:,t], agrid, bgrid, zgrid)

    # Bond Choice
    for i = 1:size(yit,1)
      wp_t[i, t] = wp_int[w_t[i, t], pension[i]]
      c_t[i, t] = xt - wp_t[i, t]
      w_t[i, t+1] = r*wp_t[i, t]

      c_t[i, t] < 0.0 ? negcons +=  1 : 0
    end
  end
  negcons == 0 || println("\t$negcons negative consumption choices!")

  avgy = mean(yit); w_t /= avgy; c_t /= avgy
  percentiles = get_percentiles(w_t)

  return c_t, w_t, wp_t, percentiles
end
