################################################################################
############################# GRID CONSTRUCTION ################################
################################################################################

function grids(s_f_i::Array{Float64,3}, stdy::Array{Float64,1},
  xpoints::T, apoints::T, bpoints::T, zpoints::T, wpoints_R::T, ypoints_R::T,
  power::Float64, r::Float64, g_t::Array{Float64,1}, pension::Array{Float64,1},
  const_bel::Bool) where T<:Integer

  println("Construct Grids")
  tW = size(s_f_i,3)

  ybelief = Array{Float64}(undef, size(s_f_i,2), tW)
  ymaxbelief = Array{Float64}(undef, tW, 1); yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp(g_t[t] + [1.0, t, 1]'*s_f_i[:, i, t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # WEALTH GRID #
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  xmin = Array{Float64}(undef, tW);
  xmax = 5*ymaxbelief
  xmin[tW] = -2*yminbelief[tW]
  for t = (tW-1):-1:1
    xmin[t] = xmin[t+1]/r - yminbelief[t]
  end

  xgrid = Array{Float64}(undef, xpoints, tW); xgridexp = similar(xgrid)
  for t = tW:-1:1
    xdistexp = (xmax[t] - xmin[t])^(1/power)
    xinc = xdistexp/(xpoints-1)
    for i = 1:xpoints
      xgridexp[i, t] = (i-1)*xinc
    end
    xgrid[:, t] = xgridexp[:, t].^power .+ xmin[t]
  end

  # BELIEF GRIDS #
  (std(s_f_i[1,:,:]) > 0.01) || (apoints = 3)
  (std(s_f_i[2,:,:]) > 0.01) || (bpoints = 3)

  if const_bel
    agrid = collect(range(minimum(s_f_i[1, :, 2:tW]), stop=maximum(s_f_i[1, :, :]), length=apoints))
    bgrid = collect(range(minimum(s_f_i[2, :, :]), stop=maximum(s_f_i[2, :, :]), length=bpoints))
    zgrid = collect(range(minimum(s_f_i[3, :, :]), stop=maximum(s_f_i[3, :, :]), length=zpoints))
  else
    agrid = Array{Float64}(undef, apoints,tW)
    bgrid = Array{Float64}(undef, bpoints,tW)
    zgrid = Array{Float64}(undef, zpoints,tW)
    for t = 1:tW
      agrid[:, t] = collect(range(minimum(s_f_i[1, :, t]), stop=maximum(s_f_i[1, :, t]), length=apoints))
      bgrid[:, t] = collect(range(minimum(s_f_i[2, :, t]), stop=maximum(s_f_i[2, :, t]), length=bpoints))
      zgrid[:, t] = collect(range(minimum(s_f_i[3, :, t]), stop=maximum(s_f_i[3, :, t]), length=zpoints))
    end
  end

  if std(s_f_i[1,:,:]) < 0.01
    agrid = collect(range(minimum(s_f_i[1, :, 2:tW])-0.001, stop=maximum(s_f_i[1, :, :])+0.001, length=apoints))
    bgrid = collect(range(minimum(s_f_i[2, :, :])-0.001, stop=maximum(s_f_i[2, :, :])+0.001, length=bpoints))
  end

  # RETIREMENT GRIDS #
  wminR = 1.05*xgrid[1,end]
  wmaxR = 0.75*xgrid[end,end]
  wgrid_R = Array{Float64}(undef, wpoints_R)
  wgridexp_R = Array{Float64}(undef, wpoints_R)
  wdistexp = (wmaxR - wminR)^(1/power)
  winc = wdistexp/(wpoints_R-1)
  for i = 1:wpoints_R
    wgridexp_R[i] = (i-1)*winc
  end
  wgrid_R = wgridexp_R.^power .+ wminR

  ygrid_R = collect(range(0.8*minimum(pension), stop=maximum(pension), length=ypoints_R))

  println("\tC-i-h grid: $(round.([xgrid[1,1] xgrid[end,1]],digits=2)) in period 1, "*
    "$(round.([xgrid[1,end] xgrid[end,end]],digits=2)) in period 40")
  println("\tBelief grids: α $(round.([agrid[1] agrid[end]],digits=2)), β "*
    "$(round.([bgrid[1] bgrid[end]],digits=2)), z $(round.([zgrid[1] zgrid[end]],digits=2))")
  println("\tRetirement grids: w_R $(round.([wgrid_R[1,1] wgrid_R[end,1]],digits=2)), "*
    "y_R $(round.([ygrid_R[1,1] ygrid_R[end,1]],digits=2))")

  return xgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end
