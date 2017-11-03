setwd( "./Laneri")
source("../lhs.R")
source("../priors-lhs.R")
dir_model = "."
#dir.create(file.path(dir_model,"mcmc"))
#dir.create(file.path(dir_model,"lhs"))
dir_lhs=file.path(dir_model,"lhs")
dir_bin=file.path(dir_model,"bin")
dir_det=file.path(dir_model,"mcmc")


# BUILD MODEL
#cmd <- sprintf("cd %s ; ssm", dir_model)
#system(cmd)

# PERFORM SIMPLEX ON LHS SAMPLE
#generate_lhs(10, priors_laneri)
#do_lhs(dir_lhs,dir_bin)
#summarize_lhs(dir_lhs,dir_model) # to be executed alone


# PERFORM ADAPTIVE MCMC
#cmd <- sprintf("cd %s/bin; cat ../theta_map_simplex.json | ./pmcmc -J 1 -M 100000 --trace --traj --hat -I 1 --root ../%s > ../%s/theta_pmcmc_det.json",dir_model,dir_det, dir_model)
#system(cmd)
#cmd <- sprintf("cd %s/bin; cat ../theta_pmcmc_det.json | ./pmcmc -J 1 -M 100000 --trace --traj --hat -I 1 --root ../%s > ../%s/theta_pmcmc_det2.json",dir_model,dir_det, dir_model)
#system(cmd)
#cmd <- sprintf("cd %s/bin; cat ../theta_pmcmc_det2.json | ./pmcmc -J 1 -M 100000  --trace --traj --hat -I 1 --root ../%s > ../%s/theta_pmcmc_det3.json",dir_model,dir_det, dir_model)
#system(cmd)

# PERFORM MCMC (without adaptation)
#cmd <- sprintf("cd %s/bin; cat ../theta_pmcmc_det3.json | ./pmcmc -J 1 -M 100000 --switch 500000   --trace --traj --hat -I 1 --root ../%s > ../%s/theta_pmcmc_det4.json",dir_model,dir_det, dir_model)
#system(cmd)

# SIMULATE FROM POSTERIOR
#cmd <- sprintf(" cd %s/bin; ssm-predict ../theta_pmcmc_det4.json ../mcmc/X_1.csv ../mcmc/trace_1.csv 2002-01-07 | ./simul -I 128 --start 2002-01-07  --end 2016-01-01 --verbose --freq 7 --traj --root ../%s",dir_model,dir_det)
#system(cmd)