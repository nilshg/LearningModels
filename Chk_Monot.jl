# Check monotonicity

function checkmonotonicity(v::Array{Float64, 6}, t::Int64)

  vnow = v[:, :, :, :, :, t]
  vsize = size(vnow)
  dim1count = 0
  dim2count = 0
  dim3count = 0
  dim4count = 0
  dim5count = 0
  # Monotonicity in first dimension
  for f = 1:vsize[5]
    for d = 1:vsize[4]
      for c = 1:vsize[3]
        for b = 1:vsize[2]
          dim1count += ~issorted(vnow[:, b, c, d, f])
        end
      end
    end
  end
  # Monotonicity in second dimension
  for f = 1:vsize[5]
    for d = 1:vsize[4]
      for c = 1:vsize[3]
        for a = 1:vsize[1]
          dim2count += issorted(vnow[a, :, c, d, f])
        end
      end
    end
  end
  # Monotonicity in third dimension
  for f = 1:vsize[5]
    for d = 1:vsize[4]
      for b = 1:vsize[2]
        for a = 1:vsize[1]
          dim3count += ~issorted(vnow[a, b, :, d, f])
        end
      end
    end
  end
  # Monotonicity in fourth dimension
  for f = 1:vsize[5]
    for c = 1:vsize[3]
      for b = 1:vsize[2]
        for a = 1:vsize[1]
          dim4count += ~issorted(vnow[a, b, c, :, f])
        end
      end
    end
  end
  # Monotonicity in fifth dimension
  for d = 1:vsize[4]
    for c = 1:vsize[3]
      for b = 1:vsize[2]
        for a = 1:vsize[1]
          dim5count += ~issorted(vnow[a, b, c, d, :])
        end
      end
    end
  end

  return dim1count, dim2count, dim3count, dim4count, dim5count
end


function checkmonotonicity(v::Array{Float64, 5}, t::Int64)

  vnow = v[:, :, :, :, t]
  vsize = size(vnow)
  dim1count = 0
  dim2count = 0
  dim3count = 0
  dim4count = 0
  # Monotonicity in first dimension
  for d = 1:vsize[4]
    for c = 1:vsize[3]
      for b = 1:vsize[2]
        dim1count += ~issorted(vnow[:, b, c, d])
      end
    end
  end
  # Monotonicity in second dimension
  for d = 1:vsize[4]
    for c = 1:vsize[3]
      for a = 1:vsize[1]
        dim2count += ~issorted(vnow[a, :, c, d])
      end
    end
  end
  # Monotonicity in third dimension
  for d = 1:vsize[4]
    for b = 1:vsize[2]
      for a = 1:vsize[1]
        dim3count += ~issorted(vnow[a, b, :, d])
      end
    end
  end
  # Monotonicity in fourth dimension
  for c = 1:vsize[3]
    for b = 1:vsize[2]
      for a = 1:vsize[1]
        dim4count += ~issorted(vnow[a, b, c, :])
      end
    end
  end

  return dim1count, dim2count, dim3count, dim4count
end



function checkmonotonicity(v::Array{Float64, 4}, t::Int64)

  vnow = v[:, :, :, t]
  vsize = size(vnow)
  dim1count = 0
  dim2count = 0
  dim3count = 0
  # Monotonicity in first dimension
  for c = 1:vsize[3]
    for b = 1:vsize[2]
      dim1count += ~issorted(vnow[:, b, c])
    end
  end
  # Monotonicity in second dimension
  for c = 1:vsize[3]
    for a = 1:vsize[1]
      dim2count += ~issorted(-vnow[a, :, c])
    end
  end
  # Monotonicity in third dimension
  for b = 1:vsize[2]
    for a = 1:vsize[1]
      dim3count += ~issorted(vnow[a, b, :])
    end
  end

  return dim1count, dim2count, dim3count
end


function checkmonotonicity(v::Array{Float64, 3}, t::Int64)

  vsize = size(v)
  dim1count = 0
  dim2count = 0
  # Monotonicity in first dimension
  for b = 1:vsize[2]
    dim1count += ~issorted(v[:, b, t])
  end

  # Monotonicity in second dimension
  for a = 1:vsize[1]
    dim2count += ~issorted(v[a, :, t])
  end

  return dim1count, dim2count
end
