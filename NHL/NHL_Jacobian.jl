using StatsBase

function jacobian{T:<AbstractFloat}(dopt::T, sopt::T)

  function u(c::Float64, σ=sopt)
    if c > 0.00001
       ut = c^(1-σ)/(1-σ)
     else
       ut = -10000. - 100*c^2.
     end
  end

  v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, dopt, sopt, tR)
  v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid,
                                      bgrid, zgrid, yit, g_t, r, dopt)
  v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                          stdy, k, r, dopt, ρ, c_over_x, g_t)
  c_t, w_t, wp_t = sim(wp, xgrid, agrid, bgrid,  zgrid, wgrid_R, yit, s_f_i,
                                                    pension, r, dopt, sopt, tR)

  age_pctl = [mean(w_t[:,6:15],2); mean(w_t[:,16:25],2); mean(w_t[:,26:35],2)]
  mopt = zeros(81,3)
  for i = 1:81
    mopt[i,1] = percentile(age_pctl[:,1], i)
    mopt[i,2] = percentile(age_pctl[:,2], i)
    mopt[i,3] = percentile(age_pctl[:,3], i)
  end
  mopt = collect(mopt)

  perturb = [dopt+0.005 sopt;
             dopt-0.005 sopt;
             dopt sopt+0.05;
             dopt sopt-0.05]
  pertresults = zeros(length(mopt),4)

  for i = 1:4
    function u(c::Float64, σ=sopt)
      if c > 0.00001
         ut = c^(1-σ)/(1-σ)
       else
         ut = -10000. - 100*c^2.
       end
    end

    v_R, wp_R = solveRetirement(wgrid_R, ygrid_R, r, dopt, sopt, tR)
    v, wp, c_over_x = solveTransition(v_R, wgrid_R, ygrid_R, xgrid, agrid,
                                        bgrid, zgrid, yit, g_t, r, dopt)
    v, wp, c_over_x = solveWorkingLife(v, wp, xgrid, agrid, bgrid, zgrid,
                                            stdy, k, r, dopt, ρ, c_over_x, g_t)
    c_t, w_t, wp_t = sim(wp, xgrid, agrid, bgrid,  zgrid, wgrid_R, yit, s_f_i,
                                                      pension, r, dopt, sopt, tR)

    age_pctl = [mean(w_t[:,6:15],2); mean(w_t[:,16:25],2); mean(w_t[:,26:35],2)]
    pertopt = zeros(81,3)
    for i = 1:81
      pertopt[i,1] = percentile(age_pctl[:,1], i)
      pertopt[i,2] = percentile(age_pctl[:,2], i)
      pertopt[i,3] = percentile(age_pctl[:,3], i)
    end
    pertresults[:,i] = collect(pertopt)
  end

  J = zeros(length(mopt),2)
  J[:,1] = 0.5*(pertresults[:,1] - mopt)/0.005
         + 0.5*(pertresults[:,2] - mopt)/(-0.005)
  J[:,2] = 0.5*(pertresults[:,3] - mopt)/0.05
         + 0.5*(pertresults[:,4] - mopt)/(-0.05)

  stderr = sqrt(diag(J'*J))
end
