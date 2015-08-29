###################################################################################
#########################    RETIREMENT PROBLEM      ##############################
###################################################################################

function solveRetirement{T<:AbstractFloat}(wgrid_R::Array{T,2},
  hgrid_R::Array{T,2}, ygrid_R::Array{T,1}, r::T, δ::T, λ::T)

  @printf "4. Solving the retirement problem\n"
  tR = size(wgrid_R,2)
  v_R = SharedArray(Float64,
    (size(wgrid_R,1), size(hgrid_R,1), size(ygrid_R,1), tR), pids = procs())
  wp_R = similar(v_R)

  # Value of last period of retirement
  for w = 1:size(wgrid_R,1), h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
    wt = wgrid_R[w, tR]; ht = hgrid_R[h, tR]
    yt = ygrid_R[y]
    xt = wt + yt
    v_R[w, h, y, tR] = u_h(xt, ht)
    wp_R[w, h, y, tR] = 0.0
  end

  # Compute the period tR-1 solution exactly:
  wmin = wgrid_R[1, tR]
  for w = 1:size(wgrid_R,1), h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
    wt = wgrid_R[w, tR-1]; ht = hgrid_R[h, tR-1]
    yt = ygrid_R[y]
    xt = wt + yt

    Blmn(wp::Float64) = -(u_h(xt - wp/r, ht)
                          + δ*u_h(wp+yt, (1-λ)*h + λ*(xt-wp/r)))

    optimum = optimize(Blmn, wmin, xt-0.01)

    v_R[w, h, y, tR-1] = -(optimum.f_minimum)
    wp_R[w, h, y, tR-1] = optimum.minimum
  end

  mv = checkmonotonicity(sdata(v_R[:,:,:,tR-1]))
  println("\tMonotonicity violations: [w,h,y],t=$(mv), $(tR-1)")

  for t = (tR-2):-1:1
    # INTERPOLATION
    v_R_interpol =
      interpolateV(sdata(v_R[:,:,:,t+1]), wgrid_R[:,t], hgrid_R[:,t], ygrid_R)

    # MAXIMIZATION
    wmin = wgrid_R[1, t+1]
    @inbounds @sync @parallel for w = 1:size(wgrid_R,1)
      for h = 1:size(hgrid_R,1), y = 1:size(ygrid_R,1)
        wt = wgrid_R[w, t]; ht = hgrid_R[h, t]; yt = ygrid_R[y]

        wt + yt > wmin/r ? 0 : error("xt < wmin/r!")

        (wp_R[w, h, y, t], v_R[w, h, y, t]) =
          bellOpt_R(wt, ht, yt, wmin, v_R_interpol, r, δ, λ)

        isnan(v_R[w, h, y, t]) ? nan_counter += 1 : 0
        wp_R[w, h, y, t] > wt + yt ? neg_c_counter += 1 : 0
      end
    end
    mv = checkmonotonicity(sdata(v_R[:,:,:,t]))
    println("\tMonotonicity violations: [w,h,y],t=$(mv), $t")
  end

  return sdata(v_R), sdata(wp_R)
end
