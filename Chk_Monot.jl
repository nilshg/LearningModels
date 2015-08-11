# Check monotonicity

function checkmonotonicity(v::Array{Float64,5})
  counter = zeros(5)
  # Monotonicity in first dimension
  for f = 1:size(v,5), d = 1:size(v,4), c = 1:size(v,3), b = 1:size(v,2)
    counter[1] += ~issorted(v[:, b, c, d, f])
  end
  # Monotonicity in second dimension
  for f = 1:size(v,5), d = 1:size(v,4), c = 1:size(v,3), a = 1:size(v,1)
    counter[2] += ~issorted(v[a, :, c, d, f], rev=true)
  end
  # Monotonicity in third dimension
  for f = 1:size(v,5), d = 1:size(v,4), b = 1:size(v,2), a = 1:size(v,1)
    counter[3] += ~issorted(v[a, b, :, d, f])
  end
  # Monotonicity in fourth dimension
  for f = 1:size(v,5), c = 1:size(v,3), b = 1:size(v,2), a = 1:size(v,1)
    counter[4] += ~issorted(v[a, b, c, :, f])
  end
  # Monotonicity in fifth dimension
  for d = 1:size(v,4), c = 1:size(v,3), b = 1:size(v,2), a = 1:size(v,1)
    counter[5] += ~issorted(v[a, b, c, d, :])
  end

  return counter
end


function checkmonotonicity(v::Array{Float64,4})
  counter = zeros(4)
  # Monotonicity in first dimension
  for d = 1:size(v,4), c = 1:size(v,3), b = 1:size(v,2)
    counter[1] += ~issorted(v[:, b, c, d])
  end
  # Monotonicity in second dimension
  for d = 1:size(v,4), c = 1:size(v,3), a = 1:size(v,1)
    counter[2] += ~issorted(v[a, :, c, d])
  end
  # Monotonicity in third dimension
  for d = 1:size(v,4), b = 1:size(v,2), a = 1:size(v,1)
    counter[3] += ~issorted(v[a, b, :, d])
  end
  # Monotonicity in fourth dimension
  for c = 1:size(v,3), b = 1:size(v,2), a = 1:size(v,1)
    counter[4] += ~issorted(v[a, b, c, :])
  end

  return counter
end


function checkmonotonicity(v::Array{Float64,3})
  counter = zeros(3)
  # Monotonicity in first dimension
  for c = 1:size(v,3), b = 1:size(v,2)
    counter[1] += ~issorted(v[:, b, c])
  end
  # Monotonicity in second dimension
  for c = 1:size(v,3), a = 1:size(v,1)
    counter[2] += ~issorted(v[a, :, c], rev=true)
  end
  # Monotonicity in third dimension
  for b = 1:size(v,2), a = 1:size(v,1)
    counter[3] += ~issorted(v[a, b, :])
  end

  return counter
end


function checkmonotonicity(v::Array{Float64,2})
  counter = zeros(2)
  # Monotonicity in first dimension
  for b = 1:size(v,2)
    counter[1] += ~issorted(v[:, b])
  end
  # Monotonicity in second dimension
  for a = 1:size(v,1)
    counter[2] += ~issorted(v[a, :])
  end

  return counter
end
