using StatsBase

function jacobian{T<:AbstractFloat}(dopt::T, sopt::T, uopt::T, wgrid_R=wgrid_R, ygrid_R=ygrid_R,
  r=r, tR=tR, xgrid=xgrid, agrid=agrid, bgrid=bgrid, zgrid=zgrid, yit=yit, g_t=g_t,
  stdy=stdy, k=k, ρ=ρ, s_f_i=s_f_i, pension=pension, power=power)

  v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ, power)
  v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid,
                                      bgrid, zgrid, yit, g_t, r, δ, σ)
  v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                    stdy, k, r, δ, ρ, c_over_x, g_t, σ, 0.0)
  c_t, w_t, wp_t, pct = sim(wp, wp_R, xgrid, agrid, bgrid, zgrid, wgrid_R,
                                  ygrid_R, yit, s_f_i, pension, r, δ, σ, tR)

  mopt = collect(pct[10:90,:])
  perturb = [dopt+0.005 sopt uopt;
             dopt-0.005 sopt uopt;
             dopt sopt+0.05 uopt;
             dopt sopt-0.05 uopt;
             dopt sopt uopt+0.05;
             dopt sopt uopt-0.05]
  pertresults = zeros(length(mopt),6)

  for i = 1:6
    δ = perturb[i,1]; σ = perturb[i,2]; υ = perturb[i,3]
    v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, δ, σ, tR, υ)
    v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid,
                                        bgrid, zgrid, yit, g_t, r, δ, σ)
    v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                      stdy, k, r, δ, ρ, c_over_x, g_t, σ, 0.0)
    c_t, w_t, wp_t, pct = sim(wp, wp_R, xgrid, agrid, bgrid, zgrid, wgrid_R,
                                    ygrid_R, yit, s_f_i, pension, r, δ, σ, tR)

    pertresults[:,i] = collect(pct[10:90,:])
  end

  J = zeros(length(mopt),3)
  J[:,1] = 0.5*(pertresults[:,1] - mopt)/0.005
         + 0.5*(pertresults[:,2] - mopt)/(-0.005)
  J[:,2] = 0.5*(pertresults[:,3] - mopt)/0.05
         + 0.5*(pertresults[:,4] - mopt)/(-0.05)
  J[:,3] = 0.5*(pertresults[:,5] - mopt)/0.05
         + 0.5*(pertresults[:,6] - mopt)/(-0.05)

  stderr = sqrt(diag(inv(J'*J)))
end
