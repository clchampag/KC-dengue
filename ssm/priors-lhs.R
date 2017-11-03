theta_seir <- list(
  R0 =list(dist="unif",args=list(min = 0, max = 20)), 
  beta1 = list(dist="unif",args=list(min = 0, max = 1)), 
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  rep_ndss =list(dist="unif",args=list(min = 0, max = 1)), 
  phase = list(dist="unif",args=list(min = -0.5, max = 0.5)), 
  pr_I =list(dist="unif",args=list(min = 0, max = 100)),
  pr_S =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05))
)

priors_laneri <- list(
  R0 =list(dist="unif",args=list(min = 0, max = 20)), 
  beta1 =list(dist="unif",args=list(min = 0, max = 1)), 
  phase =list(dist="unif",args=list(min = -0.5, max = 0.5)),
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  rep_ndss =list(dist="unif",args=list(min = 0, max = 1)), 
  pr_S =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05)), 
  lambda0 = list(dist="unif",args=list(min = 0.000001, max = 0.001)), 
  pr_I =  list(dist="unif",args=list(min = 0, max = 100))
)


theta_pandey <- list(
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  pr_Hi =list(dist="unif",args=list(min = 0, max = 100)), 
  pr_Hs =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05)), 
  pr_Vi =list(dist="unif",args=list(min=0.000001, max=.001)), 
  R0 = list(dist="unif",args=list(min=0, max=20)),
  beta1 =list(dist="unif",args=list(min = 0.1, max = 2)), 
  rep_ndss =list(dist="unif",args=list(min = 0, max = 1)), 
  phase =list(dist="unif",args=list(min = -0.5, max = 0.5)), 
  baV = list(dist="unif",args=list(min=0.1, max=2))
)

theta_seiar <- list(
  R0 =list(dist="unif",args=list(min = 0, max = 20)), 
  beta1 = list(dist="unif",args=list(min = 0, max = 1)), 
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  rho_h =list(dist="unif",args=list(min = 0, max = 1)), 
  rho_a =list(dist="unif",args=list(min = 0, max = 1)), 
  rep_ndss = list(dist="dirac",args=list(value=1)),
  pr_S =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05)), 
  phase = list(dist="unif",args=list(min = -0.5, max = 0.5)), 
  pr_I =list(dist="unif",args=list(min = 0, max = 100))
)

theta_seir2 <- list(
  R0 =list(dist="unif",args=list(min = 0, max = 20)), 
  beta1 = list(dist="unif",args=list(min = 0, max = 1)), 
  rep_ndss = list(dist="unif",args=list(min = 0, max = 1)), 
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  pr_S =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05)), 
  pr_S1 =list(dist="unif",args=list(min = 0.1, max = 1)),
  pr_S2 =list(dist="unif",args=list(min = 0.1, max = 1)),
  phase = list(dist="unif",args=list(min = -0.5, max = 0.5)), 
  pr_I1 = list(dist="unif",args=list(min = 0, max = 100)),
  pr_I2 = list(dist="unif",args=list(min = 0, max = 100))
)

theta_seir2psi <- list(
  R0 =list(dist="unif",args=list(min = 0, max = 20)), 
  beta1 = list(dist="unif",args=list(min = 0, max = 1)), 
  rep_ndss = list(dist="unif",args=list(min = 0, max = 1)), 
  psi = list(dist="unif",args=list(min = 0, max = 3)), 
  i =list(dist="unif",args=list(min = 0, max = 10)), 
  pr_S =list(dist="truncnorm",args=list(a = 0.2, b = 1, mean=0.44, sd=0.05)), 
  pr_S1 =list(dist="unif",args=list(min = 0.1, max = 1)),
  pr_S2 =list(dist="unif",args=list(min = 0.1, max = 1)),
  phase = list(dist="unif",args=list(min = -0.5, max = 0.5)), 
  pr_I1 = list(dist="unif",args=list(min = 0, max = 100)),
  pr_I2 = list(dist="unif",args=list(min = 0, max = 100))
)

