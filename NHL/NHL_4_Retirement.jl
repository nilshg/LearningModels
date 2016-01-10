################################################################################
#########################    RETIREMENT PROBLEM      ###########################
################################################################################

function solveRetirement{T<:AbstractFloat}(wgrid_R::Array{T,1},
  ygrid_R::Array{T,1}, r::T, δ::T, σ::T, tR::Int64)

  function get_c_1{T<:AbstractFloat}(r::T, δ::T, x::T, y::T, σ::T, tR::Int64)
    numerator = 1 - 1/r*(r*δ)^(1/σ)
    denominator = 1 - (1/r*(r*δ)^(1/σ))^tR
    margprop = numerator/denominator
    pdvresources = y*((1-1/r^tR)/(1-1/r)) + x
    return margprop * pdvresources
  end

  function simul{T<:AbstractFloat}(w::T, y::T, δ::T, σ::T, r::T, tR::Int64)
    c_t = Array(Float64, tR); x_t = similar(c_t); sum_u = 0.0

    c_t[1] = get_c_1(r, δ, w, y, σ, tR)[1]
    x_t[1] = w + y

    for t = 1:tR-1
      c_t[t+1] = (r*δ)^(1/σ)*c_t[t]
      x_t[t+1] = (x_t[t] - c_t[t])*r + y
      sum_u += (δ^t)*u(c_t[t], σ)
    end
    sum_u += (δ^tR)*u(c_t[tR], σ)

    return x_t, sum_u
  end

  wp_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))
  v_R = Array(Float64, (size(wgrid_R,1), size(ygrid_R,1), 1))

  for w = 1:size(wgrid_R,1), y = 1:size(ygrid_R,1)
    wt = wgrid_R[w]; yt = ygrid_R[y]
    (wp_R[w, y, :], v_R[w, y]) = simul(wt, yt, δ, σ, r, tR)
  end
  return v_R, wp_R
end

function solveRetirement{T<:AbstractFloat}(wgrid_R::Array{T,1},
  ygrid_R::Array{T,1}, r::T, δ::T, σ::T, tR::Int64, υ::T)

  v_R = SharedArray(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))
  wp_R = SharedArray(Float64, (size(wgrid_R,1), size(ygrid_R,1), tR))

  wretgrid = Array(Float64, (length(wgrid_R), tR))
  wgridexp = Array(Float64, length(wgrid_R))
  wdistexp = wgrid_R[end]^(1/power)
  winc = wdistexp/(length(wgrid_R)-1)
  for i = 1:length(wgrid_R)
    wgridexp[i] = (i-1)*winc
  end
  wretgrid[:, tR] = wgridexp.^power
  wretgrid[:,1] = wgrid_R
  for i = 1:length(wgrid_R)
    wretgrid[i,:] = collect(linspace(wretgrid[i,1],wretgrid[i,end],tR))
  end

  @everywhere function consdec(xt::Float64, σ::Float64, r::Float64, υ::Float64)
    obj(c::Float64, xt=xt, σ=σ, r=r, υ=υ) = -(u(c,σ) + b(r*(xt-c), υ))
    opt = Optim.optimize(obj, 0.1, xt)
    return -opt.f_minimum, opt.minimum
  end

  # Value of last period of retirement
  for w = 1:size(wgrid_R,1), y = 1:size(ygrid_R,1)
    xt = wgrid_R[w] + ygrid_R[y]
    v_R[w, y, tR], wp_R[w, y, tR] = consdec(xt, σ, r, υ)
  end

  # Compute the period tR-1 solution exactly:
  @everywhere function bellOptRet_exact(xt, σ, r, ψ, wmin, υ)
    obj_p(c::Float64, xt=xt, σ=σ, r=r, υ=υ) = -(u(c,σ) + b(r*(xt-c), υ))
    Blmn(wp::Float64, xt=xt, r=r, ψ=ψ) = -( (1-ψ[tR-1])*u(xt-wp,σ)+ψ[tR-1]*b(wp,υ)
                          + δ*(-Optim.optimize(obj_p,0.1,r*wp).f_minimum))
    optimum = optimize(Blmn, wmin, xt-0.01)
    return -(optimum.f_minimum), optimum.minimum
  end

  wmin = wretgrid[1, tR]
  @inbounds @sync @parallel for w = 1:size(wgrid_R,1)
    for y = 1:size(ygrid_R,1)
      xt = wretgrid[w, tR-1] + ygrid_R[y]
      v_R[w, y, tR-1], wp_R[w, y, tR-1] = bellOptRet_exact(xt, σ, r, ψ, wmin, υ)
    end
  end

  # Remainder with Interpolation
  @everywhere function bellOptRet(vint, xt, y, σ, r, ψ, wmin, υ)
    function Blmn(wp::Float64, vint=vint, xt=xt, r=r, ψ=ψ, υ=υ)
      -( (1-ψ[tR-1])*u(xt-wp,σ)+ψ[tR-1]*b(wp, υ) + δ*vint[r*wp, y])
    end
    optimum = optimize(Blmn, wmin, xt-0.01)
    return -(optimum.f_minimum), optimum.minimum
  end

  for t = tR-2:-1:1
    vint = interpolateV(v_R[:,:,t+1], wretgrid[:,t+2], ygrid_R)

    @sync @parallel for w = 1:size(wgrid_R,1)
      for y = 1:size(ygrid_R,1)
        xt = wretgrid[w, t] + ygrid_R[y]
        v_R[w, y, t], wp_R[w, y, t] = bellOptRet(vint,xt,y, σ, r, ψ, wmin, υ)
      end
    end
  end
  return sdata(v_R), sdata(wp_R), wretgrid

end
