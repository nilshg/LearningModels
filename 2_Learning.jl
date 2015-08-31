################################################################################
############################## BELIEF FORMATION ################################
################################################################################

function learning(user::AbstractString)

  @printf "2. Import Guvenen's learning results\n"
  # Import beliefs
  path="C:/Users/"*user*"/Dropbox/QMUL/PhD/Code/Guvenen FORTRAN Code/"
  s_f_guv_org = readdlm(path*"/SNext_in.dat")
  s_f_i = Array(Float64, (3, 100000, 40))
  for t = 1:40
    for i = 1:100000
      s_f_i[:, i, t] = s_f_guv_org[i, 3*t-2:3*t]
    end
  end

  # Import P matrix
  p_f = reshape(readdlm(path*"/Pnow.dat")', 3,3,40)

  # Parameters variances, standard deviation of income
  k = Array(Float64, (3, 40))
  stdy = zeros(40)
  @inbounds for t = 1:40
    ht = [1; t; 1]
    pt = p_f[:, :, t]
    k[:, t] = pt*ht.*(ht'*pt*ht + 0.047).^(-1.0)
    stdy[t] = collect(sqrt(ht'*p_f[:, :, t]*ht + 0.047))[1]
  end

  return s_f_i, stdy, k
end

################################################################################

function learning{T<:AbstractFloat}(α::Array{T,1},β::Array{T,1},yit::Array{T,2},
  ρ::T, var_α::T, var_β::T, cov_αβ::T, var_η::T, var_ɛ::T, fpu::T)

  @printf "2. Calculate agent's beliefs\n"
  tW = size(yit,2)
  s_f_i = Array(Float64, (3, size(yit,1),tW))
  s_0_i = repmat([mean(α) ; mean(β); 0.], 1, size(yit,1))
  for i = 1:size(yit,1)
    s_0_i[2, i] = fpu*β[i] + (1-fpu)*s_0_i[2, i]
  end
  s_f_i[:,:,1] = s_0_i

  f = [1. 0. 0.; 0. 1. 0.; 0. 0. ρ]
  q = [0. 0. 0.; 0. 0. 0.; 0. 0. var_η]
  p_f = Array(Float64, (3, 3, tW))
  p_0 = [var_α     cov_αβ       0.0;    # Directly out of Guvenen's paper
         cov_αβ (1-fpu)*var_β   0.0;
          0.0       0.0        0.0885]

  # Forecast from initial beliefs
  p_f[:, :, 1] = p_0;

  # Evolution of Var-Cov-Matrix
  stdy = Array(Float64, tW); k = Array(Float64, (3, tW))
  for t = 1:tW
    ht = [1; t; 1]
    pt = p_f[:, :, t]
    k[:, t] = pt*ht.*(ht'*pt*ht + var_ɛ).^(-1.0)
    stdy[t] = [sqrt(ht'*p_f[:, :, t]*ht + var_ɛ)][1]
    if t < tW
      p_f[:, :, t+1] = f*(pt-pt*ht.*(ht'*pt*ht+var_ɛ).^(-1.0)*ht'*pt)*f' + q
      for i = 1:size(yit,1)
        s_f_i[:, i, t+1] = f*(s_f_i[:, i, t]
                            + k[:,t].*(log(yit[i, t]) - ht'*s_f_i[:, i, t]))
      end
    end
  end

  return s_f_i, stdy, k
end
