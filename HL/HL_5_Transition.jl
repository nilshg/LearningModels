###################################################################################
#########################    TRANSITION PROBLEM      ##############################
###################################################################################

function get_pension{T<:AbstractFloat}(y::T, k_0::T, k_1::T, avgy::T)
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

function solveTransition{T<:AbstractFloat}(v_R::Array{T,4}, wgrid_R::Array{T,2},
  hgrid_R::Array{T,2}, ygrid_R::Array{T,1}, xgrid::Array{T,2}, hgrid::Array{T,2},
  agrid::Array{T,1}, bgrid::Array{T,1}, zgrid::Array{T,1}, r::T, δ::T, λ::T)

  tW = size(xgrid,2)
  wp = Array(Float64, size(xgrid,1), size(hgrid,1), size(agrid,1),
             size(bgrid,1), size(zgrid,1), tW)
  v = similar(wp); c_over_x = similar(wp)

  # Predicting pension from last period income
  ybari = mean(yit, 2)[:]
  (k_0, k_1) = linreg(yit[:, tW], ybari)
  avgy = mean(yit)

  # INTERPOLATION
  valueRETIRE = interpolateV(v_R[:,:,:,1], wgrid_R[:,1], hgrid_R[:,1], ygrid_R)

  # MAXIMIZATION
  wmin = wgrid_R[1, 1]
  for a = 1:size(agrid,1), b = 1:size(bgrid,1), z = 1:size(zgrid,1)
    at = agrid[a]; bt = bgrid[a]; zt = zgrid[z]
    yt = exp(at + bt*tW + zt)
    pension = get_pension(yt,k_0,k_1,avgy)

    for x = 1:size(xgrid,1), h = 1:size(hgrid,1)
      xt = xgrid[x, tW]; ht = hgrid[h, tW]

      (wp[x, h, a, b, z, tW], v[x, h, a, b, z, tW]) =
        bellOpt_TRANS(xt, ht, pension, wmin, valueRETIRE, r, δ, λ)

      isnan(v[x,h,a,b,z,tW]) ? "V is NaN at [w,h,y]=$([wt,ht,yt])" : 0
    end
  end
  mv = checkmonotonicity(v[:,:,:,:,:,end])
  println("\tMonotonicity violations: [w,h,a,b,z]=$(mv)")

  return v, wp, c_over_x
end
