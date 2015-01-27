################################################################################
############################## BELIEF FORMATION ################################
################################################################################

function learning(α::Array, β::Array, yit::Array, ρ::Float64, var_η::Float64,
                  var_ɛ::Float64, guvenen_distribution::Bool)

  @printf "2. Calculate agent's beliefs\n"
  tW = size(yit,2)
  f = [1 0 0; 0 1 0; 0 0 ρ]
  q = [0 0 0; 0 0 0; 0 0 var_η]
  s_f_i = Array(Float64, (3, size(yit,1),tW))

  if guvenen_distribution
      s_0_i = repmat([2.0 ; mean(β); 0], 1, size(yit,1))
  else
      s_0_i = repmat([mean(α) ; mean(β); 0], 1, size(yit,1))
  end

  p_f = Array(Float64, (3, 3, tW))
  p_0 = [0.005 -0.0002 0; -0.0002 0.0001 0; 0 0 0.0885]  # Directly out of Guvenen's paper

  # Add some prior knowledge on beta
  for i = 1:size(yit,1)
      s_0_i[2, i] = 0.65*β[i] + 0.35*s_0_i[2, i];
  end

  # Forecast from initial beliefs
  s_f_i[:,:,1] = s_0_i
  p_f[:, :, 1] = p_0;

  # Evolution of Var-Cov-Matrix
  @inbounds for t = 1:tW-1
      ht = [1; t; 1]
      pt = p_f[:, :, t]
      p_f[:, :, t+1] = f*(pt-pt*ht.*(ht'*pt*ht+var_ɛ).^(-1.0)*ht'*pt)*f' + q
  end

  # Calculate Standard Deviation (needed later on)
  stdy = Array(Float64, tW)
  for t = 1:tW
    ht = [1; t; 1]
    stdy[t] = [sqrt(ht'*p_f[:, :, t]*ht + var_ɛ)][1]
  end

  # Calculate Beliefs
  @inbounds for t = 1:tW-1
    pt = p_f[:, :, t]
    ht = [1; t; 1]
    for i = 1:size(yit,1)
      s_f_i[:, i, t+1] = f*(s_f_i[:, i, t] + pt*ht.*(ht'*pt*ht + var_ɛ).^(-1.0)
                            .*(log(yit[i, t]) - ht'*s_f_i[:, i, t]) )
    end
  end

  return s_f_i, stdy
end
