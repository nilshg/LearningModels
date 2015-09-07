################################################################################
#########################    RETIREMENT PROBLEM      ###########################
################################################################################

function solveRetirement{T<:AbstractFloat}(wgrid_R::Array{T,2},
  ygrid_R::Array{T,1}, r::T, δ::T, σ::T)

  @printf "4. Solving the retirement problem\n"
  tR = size(wgrid_R,2)
  function get_c_1{T<:AbstractFloat}(r::T, δ::T, x::T, y::T, σ::T, tR::Int64)
    numerator = 1 - 1/r*(r*δ)^(1/σ)
    denominator = 1 - (1/r*(r*δ)^(1/σ))^tR
    margprop = numerator/denominator
    pdvresources = y*((1-1/r^tR)/(1-1/r)) + x
    return margprop * pdvresources
  end

  function simul{T<:AbstractFloat}(x::T, y::T, δ::T, σ::T, r::T, tR::Int64)
    c_t = Array(Float64, tR); x_t = similar(c_t); sum_u = 0.0

    c_t[1] = get_c_1(r, δ, x, y, σ, tR)[1]
    x_t[1] = x + y

    for t = 1:tR-1
      c_t[t+1] = (r*δ)^(1/σ)*c_t[t]
      x_t[t+1] = (x_t[t] - c_t[t])*r + y
      sum_u += (δ^t)*u(c_t[t])
    end
    sum_u += (δ^tR)*u(c_t[tR])

    return x_t, sum_u
  end

  wp_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))
  v_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), 1))

  for w = 1:size(wgrid_R,1), y = 1:size(ygrid_R,1)
    xt = wgrid_R[w]; yt = ygrid_R[y]
    (wp_R[w, y, :], v_R[w, y]) = simul(xt, yt, δ, σ, r, tR)
  end
  return v_R, wp_R
end
