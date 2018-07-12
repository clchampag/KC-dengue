require(data.table)
library(ggplot2)
library(gridExtra)  
library(rjson)

dir_fig = "."
setwd("/ssm")
ROOT = getwd() 

pathSEIR = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")
datapath = paste(ROOT,"/../analysis/", sep="")


file_seir="theta_pmcmc_det2.json"
file_laneri="theta_pmcmc_det3.json"
file_pandey="theta_pmcmc_det3.json"
file_seiar="theta_pmcmc_det3.json"
file_seir2="theta_pmcmc_det3.json"
file_seir2psi="theta_pmcmc_det3.json"

# GET INDICATORS FROM JSON FILE
get_dic=function(path, file_name){
  DIC=round(fromJSON(file=file.path(path,file_name))$resources[3][[1]]$data$DIC)
}
get_nbparam=function(path, file_name){
  n_parameters=round(fromJSON(file=file.path(path,file_name))$resources[3][[1]]$data$n_parameters)
}
get_ndata=function(path, file_name){
  n_parameters=round(fromJSON(file=file.path(path,file_name))$resources[3][[1]]$data$n_data)
}

nbp_seir=get_nbparam(pathSEIR,file_seir)
nbp_laneri=get_nbparam(pathLaneri,file_laneri)
nbp_pandey=get_nbparam(pathPandey,file_pandey)
nbp_seiar=get_nbparam(pathSEIAR,file_seiar)
nbp_seir2=get_nbparam(pathSEIR2,file_seir2)
nbp_seir2psi=get_nbparam(pathSEIR2_psi,file_seir2psi)

nd_seir=get_ndata(pathSEIR,file_seir)
nd_laneri=get_ndata(pathLaneri,file_laneri)
nd_pandey=get_ndata(pathPandey,file_pandey)
nd_seiar=get_ndata(pathSEIAR,file_seiar)
nd_seir2=get_ndata(pathSEIR2,file_seir2)
nd_seir2psi=get_ndata(pathSEIR2_psi,file_seir2psi)

dic_seir=get_dic(pathSEIR,file_seir)
dic_laneri=get_dic(pathLaneri,file_laneri)
dic_pandey=get_dic(pathPandey,file_pandey)
dic_seiar=get_dic(pathSEIAR,file_seiar)
dic_seir2=get_dic(pathSEIR2,file_seir2)
dic_seir2psi=get_dic(pathSEIR2_psi,file_seir2psi)

#RMSE on NDSS data
mse_ndss=function(path, mcmc, start_date, end_date, cut_date, burn=0 ,pattern="X_128.csv", datapath=NULL){
  fileNames = list.files(file.path(path,mcmc), pattern = pattern, full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date","ran_incidence_ndss_obs" ,"incidence_ndss_obs","index"))
  
  fileData=list.files(file.path(path,"data"), pattern = "data.csv", full.names = T)
  if (is.null(datapath)==FALSE){fileData=list.files(datapath, pattern = "data_proj.csv", full.names = T)}
  data <- fread(fileData, sep = ",", select = c("date", "incidence_ndss_obs"))
  names(data)=c("date", "data")
  fileDf=fileDf[as.Date(fileDf$date)>=start_date & as.Date(fileDf$date)<=end_date & fileDf$index>=burn,]
  test=merge(fileDf, data, by="date")
  test$SE=(test$incidence_ndss_obs-test$data)^2
  test$SE_ran=(test$ran_incidence_ndss_obs-test$data)^2
  test$AE_ran=abs(test$ran_incidence_ndss_obs-test$data)
  
  mse=test[, mean(SE), by = c("index")]
  mse_ran=test[, mean(SE_ran), by = c("index")]
  mse_before=test[as.Date(test$date)< cut_date, mean(SE_ran), by = c("index")]
  mse_after=test[as.Date(test$date)>= cut_date, mean(SE_ran), by = c("index")]
  mae_ran=test[, mean(AE_ran), by = c("index")]
  
  print(c("with obs. noise",mean(sqrt(mse_ran$V1))))
  print(c("with obs. noise quantiles ",quantile(sqrt(mse_ran$V1), probs = c(0.025,0.5,0.975))))
  print(c("MAE with obs. noise quantiles ",quantile((mae_ran$V1), probs = c(0.025,0.5,0.975))))
  print(c("with obs. noise before ",mean(sqrt(mse_before$V1))))
  print(c("with obs. noise after",mean(sqrt(mse_after$V1))))
  return(c(mse=mean(sqrt(mse_ran$V1)),mse_q025=quantile(sqrt(mse_ran$V1), probs = c(0.025)),mse_q975=quantile(sqrt(mse_ran$V1), probs = c(0.975)), 
           mse_before=mean(sqrt(mse_before$V1)),mse_before_q025=quantile(sqrt(mse_before$V1), probs = c(0.025)),mse_before_q975=quantile(sqrt(mse_before$V1), probs = c(0.975)), 
           mse_after=mean(sqrt(mse_after$V1)),mse_after_q025=quantile(sqrt(mse_after$V1), probs = c(0.025)),mse_after_q975=quantile(sqrt(mse_after$V1), probs = c(0.975))))
}



SEIR_ndss=mse_ndss(pathSEIR,mcmc = "mcmc",cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
Laneri_ndss=mse_ndss(pathLaneri,mcmc = "mcmc", cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
Pandey_ndss=mse_ndss(pathPandey,mcmc = "mcmc", cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIAR_ndss=mse_ndss(pathSEIAR,mcmc = "mcmc", cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIR2_ndss=mse_ndss(pathSEIR2,mcmc = "mcmc", cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIR2_psi_ndss=mse_ndss(pathSEIR2_psi,mcmc = "mcmc", cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")


SEIR_proj=mse_ndss(pathSEIR,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
Laneri_proj=mse_ndss(pathLaneri,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
Pandey_proj=mse_ndss(pathPandey,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIAR_proj=mse_ndss(pathSEIAR,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIR2_proj=mse_ndss(pathSEIR2,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIR2_psi_proj=mse_ndss(pathSEIR2_psi,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")

################
#RMSE ON DENFREE DATA
mse_denfree=function(path,mcmc,pattern = "X_128.csv", ran, obs){
  fileNames = list.files(file.path(path,mcmc), pattern =pattern, full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date", "ran_incidence_denfree_obs","index"))
  
  fileData=list.files(file.path(path,"data"), pattern = "data.csv", full.names = T)
  data <- fread(fileData, sep = ",", select = c("date", "incidence_denfree_obs"))
  data=na.omit(data)
  names(data)=c("date", "data")
  test=merge(fileDf, data, by="date")
  test$SE=(test$ran_incidence_denfree_obs-test$data)^2
  mse=test[, mean(SE), by = c("index")]
  print(c("ran",mean(sqrt(mse$V1))))
  return(c(mse=mean(sqrt(mse$V1)),mse_q025=quantile(sqrt(mse$V1), probs = c(0.025)),mse_q975=quantile(sqrt(mse$V1), probs = c(0.975))))
}


mse_SEIR2=function(path,mcmc,pattern="X_128.csv"){
  fileNames = list.files(file.path(path,mcmc), pattern = pattern, full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date", "ran_incidence_D1_obs", "ran_incidence_D2_obs","index"))
  
  fileData=list.files(file.path(path,"data"), pattern = "data.csv", full.names = T)
  data <- fread(fileData, sep = ",", select = c("date", "incidence_D1_obs", "incidence_D2_obs"))
  data=na.omit(data)
  names(data)=c("date", "d1","d2")
  test=merge(fileDf, data, by="date")
  test$obs=test$d1+test$d2
  test$inc_denfree=test$ran_incidence_D1_obs+test$ran_incidence_D2_obs
  test$SE=(test$inc_denfree-test$obs)^2
  mse=test[, mean(SE), by = c("index")]
  print(c("ran",mean(sqrt(mse$V1))))
  return(c(mse=mean(sqrt(mse$V1)),mse_q025=quantile(sqrt(mse$V1), probs = c(0.025)),mse_q975=quantile(sqrt(mse$V1), probs = c(0.975))))
}


SEIR_denfree=mse_denfree(pathSEIR,mcmc="mcmc")
Laneri_denfree=mse_denfree(pathLaneri,mcmc="mcmc")
Pandey_denfree=mse_denfree(pathPandey,mcmc="mcmc")
SEIAR_denfree=mse_denfree(pathSEIAR,mcmc="mcmc")
SEIR2_denfree=mse_SEIR2(pathSEIR2,mcmc="mcmc")
SEIR2_psi_denfree=mse_SEIR2(pathSEIR2_psi,mcmc="mcmc")



write_mse=function(this_mse, my_digits=0, my_subset=""){
  
  return(paste0(round(this_mse[paste0("mse",my_subset)],digits = my_digits)," (",
                round(this_mse[paste0("mse",my_subset,"_q025.2.5%")],digits = my_digits),"-",
                round(this_mse[paste0("mse",my_subset,"_q975.97.5%")],digits = my_digits),")"))
}  


# create latex table
library(xtable)
nb_parameters=data.frame(SEIR=nbp_seir,Laneri=nbp_laneri, Pandey=nbp_pandey, SEIAR=nbp_seiar,
                         SEIR2=nbp_seir2, SEIR2psi=nbp_seir2psi)
nb_obs=c(nd_seir,nd_laneri,nd_pandey,nd_seiar,nd_seir2,nd_seir2psi)

DIC=c(dic_seir,dic_laneri,dic_pandey,dic_seiar,dic_seir2,dic_seir2psi)
nothing=c("","","","","")

# SQUARED ERROR
rmse_ndss=c(write_mse(SEIR_ndss),write_mse(Laneri_ndss),write_mse(Pandey_ndss),write_mse(SEIAR_ndss),write_mse(SEIR2_ndss),write_mse(SEIR2_psi_ndss))
rmse_ndss_before=c(write_mse(SEIR_ndss,my_subset = "_before"),write_mse(Laneri_ndss,my_subset = "_before"),write_mse(Pandey_ndss,my_subset = "_before"),write_mse(SEIAR_ndss,my_subset = "_before"),write_mse(SEIR2_ndss,my_subset = "_before"),write_mse(SEIR2_psi_ndss,my_subset = "_before"))
rmse_ndss_after=c(write_mse(SEIR_ndss,my_subset = "_after"),write_mse(Laneri_ndss,my_subset = "_after"),write_mse(Pandey_ndss,my_subset = "_after"),write_mse(SEIAR_ndss,my_subset = "_after"),write_mse(SEIR2_ndss,my_subset = "_after"),write_mse(SEIR2_psi_ndss,my_subset = "_after"))
#rmse_ndss_before=round(c(SEIR_ndss["mse_before"],Laneri_ndss["mse_before"],Pandey_ndss["mse_before"],SEIAR_ndss["mse_before"],SEIR2_ndss["mse_before"],SEIR2_psi_ndss["mse_before"]), digits=1)
#rmse_ndss_after=round(c(SEIR_ndss["mse_after"],Laneri_ndss["mse_after"],Pandey_ndss["mse_after"],SEIAR_ndss["mse_after"],SEIR2_ndss["mse_after"],SEIR2_psi_ndss["mse_after"]), digits=1)

rmse_denfree=c(write_mse(SEIR_denfree),write_mse(Laneri_denfree),write_mse(Pandey_denfree),write_mse(SEIAR_denfree),write_mse(SEIR2_denfree),write_mse(SEIR2_psi_denfree))
rmse_proj=c(write_mse(SEIR_proj),write_mse(Laneri_proj),write_mse(Pandey_proj),write_mse(SEIAR_proj),write_mse(SEIR2_proj),write_mse(SEIR2_psi_proj))
rmse_proj_14=c(write_mse(SEIR_proj,my_subset = "_before"),write_mse(Laneri_proj,my_subset = "_before"),write_mse(Pandey_proj,my_subset = "_before"),write_mse(SEIAR_proj,my_subset = "_before"),write_mse(SEIR2_proj,my_subset = "_before"),write_mse(SEIR2_psi_proj,my_subset = "_before"))
rmse_proj_15=c(write_mse(SEIR_proj,my_subset = "_after"),write_mse(Laneri_proj,my_subset = "_after"),write_mse(Pandey_proj,my_subset = "_after"),write_mse(SEIAR_proj,my_subset = "_after"),write_mse(SEIR2_proj,my_subset = "_after"),write_mse(SEIR2_psi_proj,my_subset = "_after"))

# values
IC_MSE=rbind(nb_parameters,nb_obs,
             nothing,DIC,rmse_ndss,rmse_ndss_before,rmse_ndss_after,rmse_denfree,
             nothing,rmse_proj,rmse_proj_14,rmse_proj_15)

# row names
IC_MSE$Model=c("nb parameters","nb observations",
               "ESTIMATION SET","DIC","RMSE NDSS","RMSE NDSS 2002","RMSE NDSS 2003-2013","RMSE DENFREE",
               "TEST SET","RMSE 2014-2015","RMSE 2014","RMSE 2015")
IC_MSE=IC_MSE[c("Model","SEIR","Laneri","Pandey","SEIAR","SEIR2","SEIR2psi")]

# export table
bold <- function(x){paste0('{\\bfseries ', x, '}')}
print(xtable(IC_MSE, type = "latex"), file = file.path(dir_fig,"stat.tex"),
      include.rownames = F,hline.after = c(-1,0,2,8,12),sanitize.colnames.function = bold)

