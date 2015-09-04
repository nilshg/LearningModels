################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids(s_f_i::Array{Float64, 3}, stdy::Array, wpoints::Int64,
  apoints::Int64, bpoints::Int64, zpoints::Int64, wpoints_R::Int64,
  ypoints_R::Int64, wmaxR::Float64, power::Float64, r::Float64, tR::Int64,
  const_bel::Bool)

  @printf "3. Construct Grids\n"
  tW = size(s_f_i,3)

  ybelief = Array(Float64, (size(s_f_i,2), tW))
  ymaxbelief = Array(Float64, (tW, 1)); yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp([1.0, t, 1]'*s_f_i[:, i, t] + 3*stdy[t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  # WEALTH GRID #
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  wmin = Array(Float64, tW); wmax = similar(wmin)

  wmin[tW] = -0.7*yminbelief[tW]
  wmax[tW] = 1.5*ymaxbelief[tW]
  for t = (tW-1):-1:1
    wmin[t] = wmin[t+1]/r - 0.7*yminbelief[t]
    wmax[t] = 2*ymaxbelief[t]
  end

  wgrid = Array(Float64, (wpoints, tW)); wgridexp = similar(wgrid)

  for t = tW:-1:1
    wdistexp = (wmax[t] - wmin[t])^(1/power)
    winc = wdistexp/(wpoints-1)
    for i = 1: wpoints
      wgridexp[i, t] = (i-1)*winc
    end
    wgrid[:, t] = wgridexp[:, t].^power + wmin[t]
  end

  # BELIEF GRIDS #
  if const_bel
    agrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[1, :, 2:tW]),
                             maximum(s_f_i[1, :, :]), apoints))
    bgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[2, :, :]),
                             maximum(s_f_i[2, :, :]), bpoints))
    zgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[3, :, :]),
                             maximum(s_f_i[3, :, :]), zpoints))
  else
    agrid = Array(Float64,(apoints,tW))
    bgrid = Array(Float64,(bpoints,tW))
    zgrid = Array(Float64,(zpoints,tW))
    for t = 1:tW
      agrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[1, :, t]),
                             maximum(s_f_i[1, :, t]), apoints))
      bgrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[2, :, t]),
                             maximum(s_f_i[2, :, t]), bpoints))
      zgrid[:, t] = convert(Array{Float64,1}, linspace(minimum(s_f_i[3, :, t]),
                             maximum(s_f_i[3, :, t]), zpoints))
    end
  end

  # Adjust borrowing constraints such that lowest belief does not have an empty
  # choice set in any period:
  for t = (size(s_f_i,3)-1):-1:1
    if const_bel
      ymin = exp(agrid[1] + t*bgrid[1] + zgrid[1])
    else
      ymin = exp(agrid[1,t] + t*bgrid[1,t] + zgrid[1,t])
    end
    tightening = wgrid[1, t] - wgrid[1, t+1]/r
    while ymin + tightening < 0.01
      tightening = wgrid[1, t] - wgrid[1, t+1]/r
      wgrid[:, t] = linspace(wgrid[1, t]+0.1, wgrid[end, t], wpoints)
    end
  end

  # RETIREMENT GRIDS #
  wminR = Array(Float64, (1, tR))
  wminR[tR] = 0.1059/(0.7)   # Directly from Guvenen's code

  for t = tR:-1:2
    wminR[t-1] = wminR[t]/r + wminR[tR] - 0.02
  end

  wminR = -0.7*wminR
  wgrid_R = Array(Float64, (wpoints_R, tR))
  wgridexp_R = Array(Float64, (wpoints_R, tR))

  for t = 1:tR
    wdistexp = (wmaxR - wminR[t])^(1/power)
    winc = wdistexp/(wpoints_R-1)
    for i = 1:wpoints_R
      wgridexp_R[i, t] = (i-1)*winc
    end
    wgrid_R[:, t] = wgridexp_R[:, t].^power + wminR[t]
  end

  yminR = max(0.2*yminbelief[tW], 0.2)
  ymaxR = min(0.05*ymaxbelief[tW], 1000)
  ygrid_R = convert(Array, linspace(yminR, ymaxR, ypoints_R))

  @printf "\tWealth grid: [%.2f %.2f] in period 1, [%.2f %.2f] in period 40\n" wgrid[1,1] wgrid[end,1] wgrid[1,end] wgrid[end,end]
  @printf "\tBelief grids: α [%.2f %.2f], β [%.2f %.2f], z [%.2f %.2f]\n" agrid[1] agrid[end] bgrid[1] bgrid[end] zgrid[1] zgrid[end]
  @printf "\tRetirement grids: w_R [%.2f %.2f], y_R [%.2f %.2f]\n" wgrid_R[1, 5] wgrid_R[end, 5] ygrid_R[1] ygrid_R[end]

  return wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end

################################################################################

function grids(s_f_i::Array{Float64, 3}, wpoints::Int64, apoints::Int64,
  bpoints::Int64, zpoints::Int64, wpoints_R::Int64, ypoints_R::Int64,
  r::Float64, user::AbstractString)

  @printf "3. Construct Grids (using Guvenen's data)\n"
  path="C:/Users/"*user*"/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/"
  # WEALTH GRID #
  wgrid_org = readdlm(path*"wealth.dat")'
  wgrid = Array(Float64, (wpoints, 40)); wgridexp = similar(wgrid)

  for t = 40:-1:1
    wdistexp = (wgrid_org[end, t] - wgrid_org[1, t])^(1./power)
    winc = wdistexp/(wpoints-1)
    for i = 1: wpoints
      wgridexp[i, t] = (i-1)*winc
    end
    wgrid[:, t] = wgridexp[:, t].^2. + wgrid_org[1, t]
  end

  # BELIEF GRIDS #
  agrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[1, :, 2:end]),
                   maximum(s_f_i[1, :, :]), apoints))
  bgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[2, :, :]),
                   maximum(s_f_i[2, :, :]), bpoints))
  zgrid = convert(Array{Float64,1}, linspace(minimum(s_f_i[3, :, :]),
                   maximum(s_f_i[3, :, :]), zpoints))

  # RETIREMENT GRIDS #
  guvgrid_R_org = readdlm(path*"wealthR.dat")'
  guvgrid_R = reshape(guvgrid_R_org, 1, 7200)
  guvgrid_R = unique(reshape(guvgrid_R_org, 12, 600), 2)
  wgrid_R = Array(Float64, (wpoints_R, 30)); wgridexp = similar(wgrid_R)

  for t = 1:30
    wdistexp = sqrt(guvgrid_R[end, t] - guvgrid_R[1, t])
    winc = wdistexp/(wpoints_R-1)
    for i = 1:wpoints_R
      wgridexp[i,t] = (i-1)*winc
    end
    wgrid_R[:, t] = wgridexp[:, t].^2 + guvgrid_R[1,t]
  end

  ygrid_R = convert(Array{Float64,1}, linspace(0.24, 1000, ypoints_R))

  @printf "\tWealth grid: [%.2f %.2f] in period 1, [%.2f %.2f] in period 40\n" wgrid[1,1] wgrid[end,1] wgrid[1,end] wgrid[end,end]
  @printf "\tBelief grids: α [%.2f %.2f], β [%.2f %.2f], z [%.2f %.2f]\n" agrid[1] agrid[end] bgrid[1] bgrid[end] zgrid[1] zgrid[end]
  @printf "\tRetirement grids: w_R [%.2f %.2f], y_R [%.2f %.2f]\n" wgrid_R[1, 5] wgrid_R[end, 5] ygrid_R[1] ygrid_R[end]

  return wgrid, agrid, bgrid, zgrid, wgrid_R, ygrid_R
end
