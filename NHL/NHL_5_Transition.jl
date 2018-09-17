################################################################################
#########################    TRANSITION PROBLEM      ###########################
################################################################################

function get_pension(y, k_0, k_1, avgy)
  ytilde = (k_0 + k_1*y)/avgy

  rratio = 0.0
  if ytilde < 0.3
    rratio = 0.9*ytilde
  elseif ytilde <= 2.0
    rratio = 0.27 + 0.32*(ytilde - 0.3)
  elseif ytilde <= 4.1
    rratio = 0.814 + 0.15*(ytilde - 2.0)
  else
    rratio = 1.129
  end
  return rratio*avgy
end


function solveTransition(v_R::Array{T,3}, wgrid_R::Array{T,1}, ygrid_R::Array{T,1},
  xgrid::Array{T,2}, agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1},
  yit::Array{T,2}, g_t::Array{T,1}, r::T, δ::T, σ::T) where T<:AbstractFloat

  tW = size(xgrid,2)
  wp = Array{Float64}(undef, size(xgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1), tW)
  v = Array{Float64}(undef, size(xgrid,1), size(agrid,1), size(bgrid,1), size(zgrid,1))

  # Predicting pension from last period income
  ybari = mean(yit, dims=2)
  linreg2(x, y) = hcat(fill!(similar(x), 1), x) \ y
  (k_0, k_1) = linreg2(yit[:, tW], ybari)
  avgy = mean(yit)

  # INTERPOLATION
  valueRETIRE = interpolateV(v_R[:, :, 1], wgrid_R[:, 1], ygrid_R)

  # MAXIMIZATION
  wmin = wgrid_R[1]
  for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
    at = agrid[a]; bt = bgrid[b]; zt = zgrid[z]
    yt = exp(g_t[tW] + at + bt*tW + zt)
    pension = get_pension(yt, k_0, k_1, avgy)

    for x = 1:size(xgrid,1)
      xt = xgrid[x, tW]

      (wp[x, a, b, z, tW], v[x, a, b, z]) =
        bellOpt_TRANS(xt, pension, wmin, valueRETIRE, r, δ, σ)
    end
  end

  return v, wp
end
