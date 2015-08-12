# Check monotonicity

function checkmonotonicity(v::Array{Float64})
  counter = Array(Int64, ndims(v))

  for i = 1:ndims(v)
    if i == 2
      counter[i] = sum(mapslices(issorted, v, i))
    else
      counter[i] = sum(~mapslices(issorted, v, i))
    end
  end

  return counter
end
