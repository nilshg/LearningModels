##############################################
############ GRID CONSTRUCTION ###############
##############################################

@printf "3. Construct Grids\n"

@printf "\t3.1 Unique Income Grid\n"

# Belief grid
maxalph = zeros(1, T);
minalph = zeros(1, T);
maxbeta = zeros(1, T);
minbeta = zeros(1, T);
maxzet = zeros(1, T);
minzet = zeros(1,T);

for t = 1:T
    maxalph[1, t] = maximum(S_f_i[1, :, t])
    maxbeta[1, t] = maximum(S_f_i[2, :, t])
    maxzet[1, t] = maximum(S_f_i[3, :, t])
    minalph[t] = minimum(S_f_i[1, :, t])
    minbeta[t] = minimum(S_f_i[2, :, t])
    minzet[t] = minimum(S_f_i[3, :, t])
end;

ALPHA = zeros(Sgridpoints-1, T)
BETA = zeros(Sgridpoints-1, T)
Z = zeros(Sgridpoints-1, T)

for t = 1:T
    ALPHA[:, t] = linspace(minalph[t], maxalph[t], Sgridpoints-1)';
    BETA[:, t] = linspace(minbeta[t], maxbeta[t], Sgridpoints-1)';
    Z[:, t] = linspace(minzet[t], maxzet[t], Sgridpoints-1)';
end

# Matrix with minimum belief triplet in all columns
MINBELIEF = zeros(3, agents*bs, T)
for t = 1:T
    MINBELIEF[1, :, t] = minalph[t]
    MINBELIEF[2, :, t] = minbeta[t]
    MINBELIEF[3, :, t] = minzet[t]
end

# Matrix that measures how far the an agent's belief is away from the minimum
DIST = S_f_i - MINBELIEF

# Matrix that measures how large one bin is
DIFF = zeros(3, agents*bs, T)
for t = 1:T
    DIFF[1, :, t] = ALPHA[2, t] - ALPHA[1, t]
    DIFF[2, :, t] = BETA[2, t] - BETA[1, t]
    DIFF[3, :, t] = Z[2, t] - Z[1, t]
end

BELIEFGRID = ceil(DIST./DIFF)
BELIEFGRID = max(BELIEFGRID, 1);

Qt = zeros(T, 1)
S_f = zeros(3, 1200, 40)
for t = 1:T
    a = unique(BELIEFGRID[:, :, t]', 1)
    S_f[:, 1:size(a, 1), t] = a'
    Qt[t] = size(a, 1)
end

for t = 1:T
    for i = 1:Qt[t]
        S_f[1, i, t] = ALPHA[S_f[1, i, t], t]
        S_f[2, i, t] = BETA[S_f[2, i, t], t]
        S_f[3, i, t] = Z[S_f[3, i, t], t]
    end
end

Ybelief = zeros(ygridpoints, 1200, T)
for t = 1:T
    for j = 1:Qt[t]
        ymean = H[:, t]'*S_f[:, j, t]
        y_low = exp(ymean - 3*stdy[t])
        y_high = exp(ymean + 3*stdy[t])
        Ybelief[:, j, t] = linspace(y_low[1], y_high[1], ygridpoints)
    end
end;

# Habit grid
# Assuming that no one can build up a habit stock larger than
# 5 times the mean wage at age T we can build the habit grid as:

@printf "\t3.2 Habit Grid\n"

hgrid = zeros(hgridpoints, T)
for t = 1:T
    hgrid[:, t] = linspace(0.1, 2*maximum(Ybelief[ygridpoints, :, t]), hgridpoints)
end;


# Wealth grid
# Maximum wealth is given by three times the highest earnings
# Minumum wealth is given by some ad hoc constraint

@printf "\t3.3 Wealth Grid\n"

power = 2
wmin = zeros(1, T)
wmax = zeros(1, T)
wgrid = zeros(wgridpoints, T)
wgridexp = zeros(wgridpoints, T)

for t = 1:T
    wmin[t] = -2*minimum(Ybelief[1, 1:Qt[t], t])
    wmax[t] = 4*maximum(Ybelief[ygridpoints, :, t])

    wdistexp = (wmax[t] - wmin[t])^(1/power)
    winc = wdistexp/(wgridpoints-1)
    for i = 1: wgridpoints
        wgridexp[i, t] = (i-1)*winc
    end
    wgrid[:, t] = wgridexp[:, t].^power + wmin[t]
end

@printf "\t3.4 Grids for α, β, z\n"

agrid = zeros(agridpoints, T)
bgrid = zeros(bgridpoints, T)
zgrid = zeros(zgridpoints, T)

for t = 1:T
  agrid[:, t] = linspace(minimum(S_f[1, 1:Qt[t], t]), maximum(S_f[1, :, t]), agridpoints)
  bgrid[:, t] = linspace(minimum(S_f[2, 1:Qt[t], t]), maximum(S_f[2, :, t]), bgridpoints)
  zgrid[:, t] = linspace(minimum(S_f[3, 1:Qt[t], t]), maximum(S_f[3, :, t]), zgridpoints)
end

@printf "\3.5 Retirement Grids\n"
wminR = zeros(1, TR)
wminR[TR] = -0.1059   # Directly from Guvenen's code

for t = TR:-1:2
  wminR[t-1] = wminR[t]/R + wminR[TR]
end

wminR[1:(TR-1)] = -0.7*wminR[1:(TR-1)]
wgrid_R = zeros(wgridpoints, TR)
wgridexp = zeros(wgridpoints, TR)
hgrid_R = zeros(hgridpoints, TR)
power = 3

for t = 1:TR
  wdistexp = (wmaxR - wminR[t])^(1/power)
  winc = wdistexp/(wgridpoints-1)
  for i = 1:wgridpoints
    wgridexp[i, t] = (i-1)*winc
  end
  wgrid_R[:, t] = wgridexp[:, t].^power + wminR[t]
  hgrid_R[:, t] = linspace(wgrid_R[1, t], 0.5*wgrid_R[end, t], hgridpoints)
end
hgrid_R = max(hgrid_R, 0.1)

yminR = max(0.2*minimum(Ybelief[:, :, T]), 0.2)
ymaxR = min(0.2*maximum(Ybelief[:, :, T]), 1000)
ygrid_R = linspace(yminR, ymaxR, ygridpoints_R)
