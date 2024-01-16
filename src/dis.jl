using Distributions
lambda_noise = 0.005
dt = 0.1

for i in 1:10000
  a=rand(Poisson(lambda_noise * dt))
  if a > 0
  println(a)
  end

end

println(Poisson(lambda_noise * dt))