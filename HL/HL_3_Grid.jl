################################################################################
##########################    GRID CONSTRUCTION    #############################
################################################################################

function Grids{T<:Int}(S_f_i::Array, stdy::Array, agents::T, bs::T, power::Float64,
               wgridpoints::T, hgridpoints::T, agridpoints::T, bgridpoints::T,
               zgridpoints::T, wgridpoints_R::T, hgridpoints_R::T,
               ygridpoints_R::T, wmaxR::Float64, R::Float64, TW::T, TR::T)

  @printf "3. Construct Grids\n"

  @printf "\t3.1 Unique Income Grid\n"
  H = [ones(1,TW); linspace(1,TW,TW)'; ones(1,TW)]
  Ybelief = Array(Float64, (agents*bs, TW))
  Ymaxbelief = Array(Float64, (TW, 1))
  Yminbelief = Array(Float64, (TW, 1))
  for t = 1:TW
    for i = 1:agents*bs
      Ybelief[i, t] = exp(H[:, t]'*S_f_i[:, i, t] + 3*stdy[t])[1]
    end
    Ymaxbelief[t] = maximum(Ybelief[:, t])*exp(3*stdy[t])
    Yminbelief[t] = minimum(Ybelief[:, t])/exp(3*stdy[t])
  end

  # Habit grid
  # Assuming that no one can build up a habit stock larger than
  # 5 times the mean wage at age TW we can build the habit grid as:

  @printf "\t3.2 Habit Grid\n"
  hgrid = Array(Float64, (hgridpoints, TW))
  for t = 1:TW
      hgrid[:, t] = linspace(0.1, 2*Ymaxbelief[t], hgridpoints)
  end

  # Wealth grid
  # Maximum wealth is given by three times the highest earnings
  # Minumum wealth is given by some ad hoc constraint

  @printf "\t3.3 Wealth Grid\n"
  wmin = Array(Float64, (1, TW))
  wmax = Array(Float64, (1, TW))
  wgrid = Array(Float64, (wgridpoints, TW))
  wgridexp = Array(Float64, (wgridpoints, TW))

  for t = 2:TW
    wmin[t] = -1.2*Yminbelief[t]
    wmax[t] = 2*Ymaxbelief[t]
    wdistexp = (wmax[t] - wmin[t])^(1/power)
    winc = wdistexp/(wgridpoints-1)
    for i = 1: wgridpoints
        wgridexp[i, t] = (i-1)*winc
    end
    wgrid[:, t] = wgridexp[:, t].^power + wmin[t]
  end

  ##############################
  #### ATTENTION: HACK!!!!!! ###
  ##############################
  wgrid[:, 1] = wgrid[:, 2]

  @printf "\t3.4 Grids for α, β, z\n"
  agrid = linspace(minimum(S_f_i[1, 2:TW, :]), maximum(S_f_i[1, :, :]), agridpoints)
  bgrid = linspace(minimum(S_f_i[2, :, :]), maximum(S_f_i[2, :, :]), bgridpoints)
  zgrid = linspace(minimum(S_f_i[3, :, :]), maximum(S_f_i[3, :, :]), zgridpoints)


  @printf "\t3.5 Retirement Grids\n"
  wminR = zeros(1, TR)
  wminR[TR] = -0.1059/7   # Directly from Guvenen's code

  for t = TR:-1:2 # This is not quite the same but gets the minimum about right
    wminR[t-1] = wminR[t]/R + 0.1*Yminbelief[TW] - 0.02
  end

  wminR = -0.7*wminR

  wgrid_R = zeros(wgridpoints_R, TR)
  wgridexp = zeros(wgridpoints_R, TR)
  hgrid_R = zeros(hgridpoints_R, TR)

  for t = 1:TR
    wdistexp = (wmaxR - wminR[t])^(1/power)
    winc = wdistexp/(wgridpoints_R-1)
    for i = 1:wgridpoints
      wgridexp[i, t] = (i-1)*winc
    end
    wgrid_R[:, t] = wgridexp[:, t].^power + wminR[t]
    hgrid_R[:, t] = linspace(wgrid_R[1, t], 0.5*wgrid_R[end, t], hgridpoints_R)
  end
  hgrid_R = max(hgrid_R, 0.1)

  yminR = max(0.2*Yminbelief[TW], 0.2)
  ymaxR = min(0.2*Ymaxbelief[TW], 1000)
  ygrid_R = linspace(yminR, ymaxR, ygridpoints_R)

  @printf "\t3.6 Check budget constraints for feasibility\n"
  for t = 1:TW-1
    wt = wgrid[1, t]
    yt = exp(agrid[1] + bgrid[1]*t + zgrid[1])
    wmin = wgrid[1, t+1]
    if wt + yt - wmin/R < 0.01
      @printf "\tError: Cash on hand is too low, at t = %d\n" t
    end
  end

  if wgrid[1, TW] + exp(agrid[1] + bgrid[1]*TW + zgrid[1]) < wgrid_R[1, 1]/R
    @printf "\tError: Cash on hand is too low in last period of working life\n"
  end

  return wgrid, hgrid, agrid, bgrid, zgrid, wgrid_R, hgrid_R, ygrid_R
end
