x = [i for i = 1:100000.0]
y = [i for i = 1:10.0]

f(x::Float64, y::Float64) = x*y

result = Array(Float64, (length(x), length(y)))

function solve()
  for y_now = 1:size(y, 1)
    for x_now = 1:size(x, 1)
      @inbounds result[x_now, y_now] =  f(x[x_now], y[y_now])
    end
  end
  return result
end
solve()
@time solve()
# Execution time: 0.28s, 127837464 bytes allocated, ~20% gc time

function solve(result::Array{Float64, 2})
  for y_now = 1:size(y, 1)
    for x_now = 1:size(x, 1)
      @inbounds result[x_now, y_now] =  f(x[x_now], y[y_now])
    end
  end
  return result
end
solve(result)
@time solve(result)
# Execution time: 0.28s, 127837464 bytes allocated, ~20% gc time

function solve(result::Array{Float64, 2}, y::Array{Float64, 1},
               x::Array{Float64, 1})
  for y_now = 1:size(y, 1)
    for x_now = 1:size(x, 1)
      @inbounds result[x_now, y_now] =  f(x[x_now], y[y_now])
    end
  end
  return result
end
solve(result, y, x)
@time solve(result, y, x)
# Execution time: 0.0009s, 80 bytes allocated

function solve(result::Array{Float64, 2}, y::Array{Float64, 1},
               x::Array{Float64, 1}, f::Function)
  for y_now = 1:size(y, 1)
    for x_now = 1:size(x, 1)
      @inbounds result[x_now, y_now] =  f(x[x_now], y[y_now])
    end
  end
  return result
end
solve(result, y, x, f)
@time solve(result, y, x, f)
# Execution time: 0.078s, 63918320 bytes allocated, ~29% gc time
