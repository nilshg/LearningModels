################################################################################
############################## BELIEF FORMATION ################################
################################################################################

function Learning(agents::Int64, bs::Int64, TW::Int64, alpha::Array, beta::Array,
                  Yit::Array, ρ::Float64, var_eta::Float64, var_eps::Float64,
                  init_z::Float64)

  @printf "2. Calculate agent's beliefs\n"

  H = [ones(1,TW); linspace(1,TW,TW)'; ones(1,TW)]
  F = [1 0 0; 0 1 0; 0 0 ρ]
  Q = [0 0 0; 0 0 0; 0 0 var_eta]
  S_f_i = zeros(3, agents*bs,TW)
  S_0_i = repmat([2.0 ; mean(beta); 0], 1, agents*bs)
  P_f = zeros(3, 3, TW)
  P_0 = [0.005 -0.0002 0; -0.0002 0.0001 0; 0 0 0.0885]  # Directly out of Guvenen's paper

  # Add some prior knowledge on beta
  for i = 1:agents*bs
      S_0_i[2, i] = 0.65*beta[i] + 0.35*S_0_i[2, i];
  end

  # Forecast from initial beliefs
  S_f_i[:,:,1] = S_0_i
  P_f[:, :, 1] = P_0;

  # Evolution of Var-Cov-Matrix
  for t = 1:TW-1
      Ht = H[:, t]
      Pf = P_f[:, :, t]
      P_f[:, :, t+1] = F*( Pf - Pf * Ht .* (Ht' * Pf * Ht + var_eps).^(-1.0) * Ht' *  Pf )*F' + Q
  end

  # Calculate Standard Deviation (needed later on)
  stdy = zeros(1,TW)
  for t = 1:TW
      std = sqrt( H[:, t]' * P_f[:, 3*t-2:3*t] * H[:, t] )
      stdy[t] = std[1]
  end

  # Calculate Beliefs
  tic()
  for t = 1:TW-1
    Pc = P_f[:, :, t]
    Hc = H[:, t]
    for i = 1:agents*bs
     S_f_i[:, i, t+1] = F*( S_f_i[:, i, t] + Pc * Hc .* (Hc' * Pc * Hc + var_eps).^(-1.0) .* (log(Yit[i, t]) - Hc' * S_f_i[:, i, t]) )
    end
  end
  @printf "\tBelief calculation took %d seconds.\n" toq()
  return S_f_i, P_f, stdy
end
