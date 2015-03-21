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
    s_f_guv_org =
      readdlm("C:\\Users\\tew207\\Dropbox\\QMUL\\PhD\\Code\\Guvenen FORTRAN code\\SNext_in.dat")

    s_f_guv = Array(Float64, (3, 100000, 40))
    for t = 1:40
      for i = 1:100000
        s_f_guv[:, i, t] = s_f_guv_org[i, 3*t-2:3*t]
      end
    end
    s_f_i[:,:,1] = s_f_guv[:, :, 1]
  else
    s_0_i = repmat([mean(α) ; mean(β); 0], 1, size(yit,1))
    for i = 1:size(yit,1)
      s_0_i[2, i] = 0.65*β[i] + 0.35*s_0_i[2, i]
    end
    s_f_i[:,:,1] = s_0_i
  end

  p_f = Array(Float64, (3, 3, tW))
  p_0 = [0.005 -0.0002 0; -0.0002 0.0001 0; 0 0 0.0885]  # Directly out of Guvenen's paper

  # Forecast from initial beliefs
  p_f[:, :, 1] = p_0;

  # Evolution of Var-Cov-Matrix
  k = Array(Float64, (3, tW))
  @inbounds for t = 1:tW-1
      ht = [1; t; 1]
      pt = p_f[:, :, t]
      k[:, t] = pt*ht.*(ht'*pt*ht + var_ɛ).^(-1.0)
    if t < tW
        p_f[:, :, t+1] = f*(pt-pt*ht.*(ht'*pt*ht+var_ɛ).^(-1.0)*ht'*pt)*f' + q
    end
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
      s_f_i[:, i, t+1] = f*(s_f_i[:, i, t]
                            + k[:,t].*(log(yit[i, t]) - ht'*s_f_i[:, i, t]))
    end
  end

  return s_f_i, stdy, k
end
