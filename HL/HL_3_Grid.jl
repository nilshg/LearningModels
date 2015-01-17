################################################################################
############################# GRID CONStRUCTION ################################
################################################################################

function grids(s_f_i::Array{Float64, 3}, stdy::Array, wpoints::Int64,
               apoints::Int64, bpoints::Int64, zpoints::Int64, wpoints_R::Int64,
               ypoints_R::Int64, wmaxR::Float64, power::Float64, r::Float64,
               tR::Int64, guv_dist::Bool, const_beliefs::Bool)

  @printf "3. Construct Grids\n"
  @printf "\t3.1 Belief Grid\n"
  tW = size(s_f_i,3)
  @printf "\t3.1 Unique Income Grid\n"
  hmat = [ones(1,tW); linspace(1,tW,tW)'; ones(1,tW)]
  ybelief = Array(Float64, (agents*bs, tW))
  ymaxbelief = Array(Float64, (tW, 1))
  yminbelief = similar(ymaxbelief)
  for t = 1:tW
    for i = 1:size(s_f_i,2)
      ybelief[i, t] = exp(hmat[:, t]'*s_f_i[:, i, t] + 3*stdy[t])[1]
    end
    ymaxbelief[t] = maximum(ybelief[:, t])*exp(3*stdy[t])
    yminbelief[t] = minimum(ybelief[:, t])/exp(3*stdy[t])
  end

  @printf "\t3.2 Habit Grid\n"
  hgrid = Array(Float64, (hpoints, tW))
  for t = 1:tW
      hgrid[:, t] = linspace(0.1, 2*ymaxbelief[t], hpoints)
  end

  # Wealth grid
  # Maximum wealth is given by three times the highest earnings
  # Minimum wealth is given by some ad hoc constraint

  @printf "\t3.2 Wealth Grid\n"
  if guv_dist
    guvgrid = readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Julia\\Guvenen\\wealth.dat")'
    wgrid = Array(Float64, (wpoints, size(s_f_i,3)))
    for t = 1:size(s_f_i,3)
      wgrid[:, t] = linspace(guvgrid[1, t], guvgrid[end, t], wpoints)
    end
  else
    wmin = Array(Float64, tW)
    wmax = Array(Float64, tW)
    wgrid = Array(Float64, (wpoints, tW))
    wgridexp = Array(Float64, (wpoints, tW))
    for t = 1:tW
      wmin[t] = -yminbelief[t]
      wmax[t] = 2*ymaxbelief[t]

      wdistexp = (wmax[t] - wmin[t])^(1/power)
      winc = wdistexp/(wpoints-1)
      for i = 1: wpoints
          wgridexp[i, t] = (i-1)*winc
      end
      wgrid[:, t] = wgridexp[:, t].^power + wmin[t]
    end
  end
  ################## HACK!! ####################
  wgrid[:, 1] = linspace(wgrid[1, 2], wgrid[end, 1], wpoints)
  ###############################################

  @printf "\t3.5 Grids for α, β, z\n"
  if const_beliefs
    agrid = linspace(minimum(s_f_i[1, :, 2:tW]), maximum(s_f_i[1, :, :]), apoints)
    bgrid = linspace(minimum(s_f_i[2, :, :]), maximum(s_f_i[2, :, :]), bpoints)
    zgrid = linspace(minimum(s_f_i[3, :, :]), maximum(s_f_i[3, :, :]), zpoints)
  else
    for t = 1:size(s_f_i,3)
      agrid[:, t] = linspace(minimum(s_f_i[1, :, t]), maximum(s_f_i[1, :, t]), apoints)
      bgrid[:, t] = linspace(minimum(s_f_i[2, :, t]), maximum(s_f_i[2, :, t]), bpoints)
      zgrid[:, t] = linspace(minimum(s_f_i[3, :, t]), maximum(s_f_i[3, :, t]), zpoints)
    end
  end

  # Adjust borrowing constraints such that lowest belief does not have an empty
  # choice set in any period:
  for t = (size(s_f_i,3)-1):-1:1
    if const_beliefs
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

  @printf "\t3.5 Retirement Grids\n"
  if guv_dist # Use Guvenen's retirement grid
    guvgrid_R_org = readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Julia\\Guvenen\\wealthR.dat")'
    guvgrid_R = reshape(guvgrid_R_org, 1, 7200)
    guvgrid_R = unique(reshape(guvgrid_R_org, 12, 600), 2)
    wgrid_R = Array(Float64, (wpoints_R, tR))
    for t = 1:30
      wgrid_R[:, t] = linspace(guvgrid_R[1, t], guvgrid_R[end, t], wpoints_R)
    end
  else
    wminR = Array(Float64, (1, tR))
    wminR[tR] = 0.1059/(0.7)   # Directly from Guvenen's code

    for t = tR:-1:2
      wminR[t-1] = wminR[t]/r + wminR[tR] - 0.02
    end

    wminR = -0.7*wminR
    wgrid_R = Array(Float64, (wpoints_R, tR))
    wgridexp_R = similar(wgrid_R)
    hgrid_R = similar(wgrid_R)

    for t = 1:tR
      wdistexp = (wmaxR - wminR[t])^(1/power)
      winc = wdistexp/(wpoints_R-1)
      for i = 1:wpoints_R
        wgridexp_R[i, t] = (i-1)*winc
      end
      wgrid_R[:, t] = wgridexp_R[:, t].^power + wminR[t]
      hgrid_R[:, t] = linspace(0.1, 0.5*wgrid_R[end, t], hpoints_R)
    end
  end

  yminR = max(0.2*yminbelief[tW], 0.2)
  ymaxR = min(0.2*ymaxbelief[tW], 1000)
  ygrid_R = linspace(yminR, ymaxR, ypoints_R)

  @printf "Wealth grid: [%.2f %.2f] in period 1, [%.2f %.2f] in period 40\n" wgrid[1,1] wgrid[end,1] wgrid[1,end] wgrid[end,end]
  @printf "Belief grids: α [%.2f %.2f], β [%.2f %.2f], z [%.2f %.2f]\n" agrid[1] agrid[end] bgrid[1] bgrid[end] zgrid[1] zgrid[end]
  @printf "Retirement grids: w_R [%.2f %.2f], y_R [%.2f %.2f]\n" wgrid_R[1, 5] wgrid_R[end, 5] ygrid_R[1] ygrid_R[end]

  return wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R ygrid_R
end
