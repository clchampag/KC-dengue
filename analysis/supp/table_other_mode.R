require(data.table)
library(ggplot2)
library(gridExtra)  
library(rjson)
library(dplyr)
dir_fig = "./supp"
setwd("/ssm")


ROOT = getwd() 
pathSEIR2 = paste(ROOT,"/SEIR2/mcmc", sep="")
pathSEIR2_paramB = paste(ROOT,"/SEIR2/mcmc-paramB", sep="")

pathSEIR2_psi = paste(ROOT,"/SEIR2-psi/mcmc", sep="")
pathSEIR2_psi_paramB = paste(ROOT,"/SEIR2-psi/mcmc-paramB", sep="")

N=161391
datapath = paste(ROOT,"/../analysis/", sep="")

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

nbp_seir2=get_nbparam(pathSEIR2,file_seir2)
nbp_seir2psi=get_nbparam(pathSEIR2_psi,file_seir2psi)
nd_seir2_other=get_ndata(pathSEIR2_paramB,file_seir2)
nd_seir2psi_other=get_ndata(pathSEIR2_psi_paramB,file_seir2psi)

nd_seir2=get_ndata(pathSEIR2,file_seir2)
nd_seir2psi=get_ndata(pathSEIR2_psi,file_seir2psi)
nd_seir2_other=get_ndata(pathSEIR2_paramB,file_seir2)
nd_seir2psi_other=get_ndata(pathSEIR2_psi_paramB,file_seir2psi)

dic_seir2=get_dic(pathSEIR2,file_seir2)
dic_seir2psi=get_dic(pathSEIR2_psi,file_seir2psi)
dic_seir2_other=get_dic(pathSEIR2_paramB,file_seir2)
dic_seir2psi_other=get_dic(pathSEIR2_psi_paramB,file_seir2psi)


#RMSE on NDSS data
mse_ndss=function(path, start_date, end_date, cut_date, burn=0 ,pattern="X_128.csv", datapath=NULL){
  fileNames = list.files(path, pattern = pattern, full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date","ran_incidence_ndss_obs" ,"incidence_ndss_obs","index"))
  
  fileData=list.files(file.path(path,"../data"), pattern = "data.csv", full.names = T)
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



SEIR2_ndss=mse_ndss(pathSEIR2, cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIR2_psi_ndss=mse_ndss(pathSEIR2_psi, cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIR2other_ndss=mse_ndss(pathSEIR2_paramB, cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")
SEIR2_psi_other_ndss=mse_ndss(pathSEIR2_psi_paramB, cut_date = "2003-01-01", start_date="2002-01-07",end_date="2013-12-30")


SEIR2_proj=mse_ndss(pathSEIR2,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01")
SEIR2_psi_proj=mse_ndss(pathSEIR2_psi,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01")
SEIR2other_proj=mse_ndss(pathSEIR2_paramB,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01")
SEIR2_psi_other_proj=mse_ndss(pathSEIR2_psi_paramB,datapath=datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01")

################
#RMSE ON DENFREE DATA

  mse_SEIR2=function(path,pattern="X_128.csv"){
    fileNames = list.files(path, pattern = pattern, full.names = T)
    fileDf = fread(fileNames[1], sep = ",", select = c("date", "ran_incidence_D1_obs", "ran_incidence_D2_obs","index"))
    
    fileData=list.files(file.path(path,"../data"), pattern = "data.csv", full.names = T)
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


SEIR2_denfree=mse_SEIR2(pathSEIR2)
SEIR2_psi_denfree=mse_SEIR2(pathSEIR2_psi)
SEIR2other_denfree=mse_SEIR2(pathSEIR2_paramB)
SEIR2_psi_other_denfree=mse_SEIR2(pathSEIR2_psi_paramB)



write_mse=function(this_mse, my_digits=0, my_subset=""){
  
  return(paste0(round(this_mse[paste0("mse",my_subset)],digits = my_digits)," (",
                round(this_mse[paste0("mse",my_subset,"_q025.2.5%")],digits = my_digits),"-",
                round(this_mse[paste0("mse",my_subset,"_q975.97.5%")],digits = my_digits),")"))
}  

############
# EPIDEMIO
hpd_trace=function(path,parameter, to_sum, my_digits=0){
  fileNames = list.files(path, pattern ="trace", full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c(parameter))
  fileDf=mutate(fileDf,param=eval(parse(text = to_sum)))
  
  my_hpd=quantile(fileDf$param, probs = c(0.5,0.025,0.975))
  my_hpd=round(my_hpd,digits = my_digits)
  return(paste0(my_hpd[1]," (",my_hpd[2],"-",my_hpd[3],")"))
}  


# Basic reproduction number
R0_S0 = function(path){
  fileNames = list.files(path, pattern = paste0("trace_1.csv"), full.names = T)
  fileTrace = fread(fileNames[1], sep = ",")
  date=seq.Date(as.Date("2002-01-14"),as.Date("2003-01-14"),by="day")
  var_R0=matrix(nrow = length(date), ncol=dim(fileTrace)[1]+1)
  var_R0[,1]=seq(1:length(date))
  for(i in (1:length(date))){
    var_R0[i,-1] = fileTrace$R0*(1 + fileTrace$beta1*sin(2*pi*(var_R0[i,1]/365+fileTrace$phase)))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_hpd_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_hpd_mean=round(my_hpd_mean,digits = 2)
  my_hpd_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_hpd_max=round(my_hpd_max,digits = 2)
  return(list(mean=paste0(my_hpd_mean[1]," (",my_hpd_mean[2],"-",my_hpd_mean[3],")"),
              max=paste0(my_hpd_max[1]," (",my_hpd_max[2],"-",my_hpd_max[3],")")))
}  


# Effective reproduction number (two strain model)
Reff_strains = function(path,N,S,S1,name){
  fileNames = list.files(path, pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",S,S1,"index"))
  names(fileDf)=c("date","S","S1","index")
  fileNames = list.files(path, pattern = paste0("trace_1.csv"), full.names = T)
  fileTrace = fread(fileNames[1], sep = ",")
  fileTrace2=fileTrace[(index*1000/length(unique(index)))%%1==0]
  fileTrace2$index=seq(0:999)-1
  
  date=data.frame(date=as.character(seq.Date(as.Date("2002-01-14"),as.Date("2013-12-31"),by="day")))
  date$x=seq(1:dim(date)[1])
  test=merge(fileDf,date)
  test=merge(test, fileTrace2, by="index")
  var_R0=matrix(nrow =  length(test$date[test$index==0]), ncol= length(unique(test$index)))
  for(i in (1:length(unique(test$x)))){
    xx = unique(test$x)[i]
    var_R0[i,] = ((test$S[test$x==xx]+test$S1[test$x==xx])/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_hpd_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_hpd_mean=round(my_hpd_mean,digits = 2)
  my_hpd_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_hpd_max=round(my_hpd_max,digits = 2)
  return(list(mean=paste0(my_hpd_mean[1]," (",my_hpd_mean[2],"-",my_hpd_mean[3],")"),
              max=paste0(my_hpd_max[1]," (",my_hpd_max[2],"-",my_hpd_max[3],")")))
}  

# Effective reproduction umber (two strain models with interaction (psi))
Reff_strains_psi = function(path,N,S,S1, I2, I12,name){
  fileNames = list.files(path, pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",S,S1,I2,I12,"index"))
  names(fileDf)=c("date","S","S1","I2","I12","index")
  fileNames = list.files(path, pattern = paste0("trace_1.csv"), full.names = T)
  fileTrace = fread(fileNames[1], sep = ",")
  fileTrace2=fileTrace[(index*1000/length(unique(index)))%%1==0]
  fileTrace2$index=seq(0:999)-1
  
  date=data.frame(date=as.character(seq.Date(as.Date("2002-01-14"),as.Date("2013-12-31"),by="day")))
  date$x=seq(1:dim(date)[1])
  test=merge(fileDf,date)
  test=merge(test, fileTrace2, by="index")
  var_R0=matrix(nrow =  length(test$date[test$index==0]), ncol= length(unique(test$index)))
  for(i in (1:length(unique(test$x)))){
    xx = unique(test$x)[i]
    var_R0[i,] =((test$S[test$x==xx]+test$psi[test$x==xx]*test$S1[test$x==xx])/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_hpd_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_hpd_mean=round(my_hpd_mean,digits = 2)
  my_hpd_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_hpd_max=round(my_hpd_max,digits = 2)
  return(list(mean=paste0(my_hpd_mean[1]," (",my_hpd_mean[2],"-",my_hpd_mean[3],")"),
              max=paste0(my_hpd_max[1]," (",my_hpd_max[2],"-",my_hpd_max[3],")")))
}  


epidemio_attack_rate=function(path, S, to_sum_S, incidence, to_sum){
  fileNames = list.files(path, pattern = "X_128.csv", full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",incidence,S,"index"))
  names(fileDf)=c("date", incidence, S,"index")
  fileDf$year=as.numeric( substr(fileDf$date, 1, 4))
  date=fileDf[, min(date), by = c("year")]
  
  names(date)=c("year","date")
  Susceptibles=merge(date,fileDf,all.x=T)
  
  Susceptibles=mutate(Susceptibles,S_tot=eval(parse(text = to_sum_S)))
  
  attack_rate=fileDf[, sum(eval(parse(text = to_sum))), by = c("year", "index")]
  tot=merge(Susceptibles,attack_rate, all=T, by=c("year","index"))
  tot$attack_rate=100*tot$V1/tot$S_tot
  tot=as.data.table(tot[tot$year<2016,])
  
  med_attack_rate=tot[, median(attack_rate), by = c("year")]
  names(med_attack_rate)=c("date","med")
  
  med=round(median(med_attack_rate$med),digits = 0)
  max=round(max(med_attack_rate$med),digits = 0)
  min=round(min(med_attack_rate$med),digits = 0)
  
  paste0(med," (",min,"-",max,")")
}

# R0
SEIR2_R0=R0_S0(pathSEIR2)
SEIR2psi_R0=R0_S0(pathSEIR2_psi)
SEIR2other_R0=R0_S0(pathSEIR2_paramB)
SEIR2psi_other_R0=R0_S0(pathSEIR2_psi_paramB)

# RE
SEIR2_D1_RE=Reff_strains(pathSEIR2,N,"S","S2",name="seir2_D1")
SEIR2_D2_RE=Reff_strains(pathSEIR2,N,"S","S1",name="seir2_D2")
SEIR2psi_D1_RE=Reff_strains_psi(path=pathSEIR2_psi,N,S="S",S1="S2",I2="I1",I12="I21",name="seir2psi_D1")
SEIR2psi_D2_RE=Reff_strains_psi(pathSEIR2_psi,N,"S","S1","I2","I12",name="seir2psi_D2")

SEIR2other_D1_RE=Reff_strains(pathSEIR2_paramB,N,"S","S2",name="seir2_D1")
SEIR2other_D2_RE=Reff_strains(pathSEIR2_paramB,N,"S","S1",name="seir2_D2")
SEIR2psi_other_D1_RE=Reff_strains_psi(path=pathSEIR2_psi_paramB,N,S="S",S1="S2",I2="I1",I12="I21",name="seir2psi_D1")
SEIR2psi_other_D2_RE=Reff_strains_psi(pathSEIR2_psi_paramB,N,"S","S1","I2","I12",name="seir2psi_D2")

# Susceptibles
SEIR2_S=hpd_trace(pathSEIR2,"pr_S",to_sum = "pr_S*100")
SEIR2psi_S=hpd_trace(pathSEIR2_psi,"pr_S",to_sum = "pr_S*100")
SEIR2other_S=hpd_trace(pathSEIR2_paramB,"pr_S",to_sum = "pr_S*100")
SEIR2psi_other_S=hpd_trace(pathSEIR2_psi_paramB,"pr_S",to_sum = "pr_S*100")

SEIR2_S1=hpd_trace(pathSEIR2,"pr_S1",to_sum = "pr_S1*100")
SEIR2psi_S1=hpd_trace(pathSEIR2_psi,"pr_S1",to_sum = "pr_S1*100")
SEIR2_S2=hpd_trace(pathSEIR2,"pr_S2",to_sum = "pr_S2*100")
SEIR2psi_S2=hpd_trace(pathSEIR2_psi,"pr_S2",to_sum = "pr_S2*100")

SEIR2other_S1=hpd_trace(pathSEIR2_paramB,"pr_S1",to_sum = "pr_S1*100")
SEIR2psi_other_S1=hpd_trace(pathSEIR2_psi_paramB,"pr_S1",to_sum = "pr_S1*100")
SEIR2other_S2=hpd_trace(pathSEIR2_paramB,"pr_S2",to_sum = "pr_S2*100")
SEIR2psi_other_S2=hpd_trace(pathSEIR2_psi_paramB,"pr_S2",to_sum = "pr_S2*100")


SEIR2_r=hpd_trace(pathSEIR2,"rep_ndss","rep_ndss*100")
SEIR2psi_r=hpd_trace(pathSEIR2_psi,"rep_ndss","rep_ndss*100")
SEIR2other_r=hpd_trace(pathSEIR2_paramB,"rep_ndss","rep_ndss*100")
SEIR2psi_other_r=hpd_trace(pathSEIR2_psi_paramB,"rep_ndss","rep_ndss*100")

psi=hpd_trace(pathSEIR2_psi,"psi","psi", my_digits = 2)
psi_other=hpd_trace(pathSEIR2_psi_paramB,"psi","psi", my_digits = 2)

AR_seir2_I=epidemio_attack_rate(pathSEIR2,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2")
AR_seir3_I=epidemio_attack_rate(pathSEIR2_psi,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2")
AR_seir2_II=epidemio_attack_rate(pathSEIR2,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I12","inc_I21"), to_sum="inc_I12+inc_I21")
AR_seir3_II=epidemio_attack_rate(pathSEIR2_psi,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I21","inc_I12"), to_sum="inc_I12+inc_I21")

AR_seir2other_I=epidemio_attack_rate(pathSEIR2_paramB,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2")
AR_seir3other_I=epidemio_attack_rate(pathSEIR2_psi_paramB,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2")
AR_seir2other_II=epidemio_attack_rate(pathSEIR2_paramB,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I12","inc_I21"), to_sum="inc_I12+inc_I21")
AR_seir3other_II=epidemio_attack_rate(pathSEIR2_psi_paramB,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I21","inc_I12"), to_sum="inc_I12+inc_I21")


############################
############################

# create latex table
library(xtable)
nb_parameters=data.frame(SEIR2_A=nbp_seir2,SEIR2_B=nbp_seir2, SEIR2psi_A=nbp_seir2psi, SEIR2psi_B=nbp_seir2psi)
nb_obs=c(nd_seir2,nd_seir2_other,nd_seir2psi,nd_seir2psi_other)

DIC=c(dic_seir2,dic_seir2_other,dic_seir2psi,dic_seir2psi_other)
nothing=c("","","","")

# SQUARED ERROR
rmse_ndss=c(write_mse(SEIR2_ndss),write_mse(SEIR2other_ndss),write_mse(SEIR2_psi_ndss),write_mse(SEIR2_psi_other_ndss))

rmse_denfree=c(write_mse(SEIR2_denfree),write_mse(SEIR2other_denfree),write_mse(SEIR2_psi_denfree),write_mse(SEIR2_psi_other_denfree))
rmse_proj=c(write_mse(SEIR2_proj),write_mse(SEIR2other_proj),write_mse(SEIR2_psi_proj),write_mse(SEIR2_psi_other_proj))

#EPIDEMIO
R0_mean=c(SEIR2_R0$mean, SEIR2other_R0$mean, SEIR2psi_R0$mean, SEIR2psi_other_R0$mean)
R0_max=c(SEIR2_R0$max, SEIR2other_R0$max, SEIR2psi_R0$max, SEIR2psi_other_R0$max)
my_psi=c("","",psi,psi_other)

HS=c(SEIR2_S,SEIR2other_S, SEIR2psi_S, SEIR2psi_other_S)
HS1=c(SEIR2_S1,SEIR2other_S1, SEIR2psi_S1, SEIR2psi_other_S1)
HS2=c(SEIR2_S2,SEIR2other_S2, SEIR2psi_S2, SEIR2psi_other_S2)

rep=c(SEIR2_r,SEIR2other_r, SEIR2psi_r, SEIR2psi_other_r)

AR=c(AR_seir2_I,AR_seir2other_I,AR_seir3_I,AR_seir3other_I)
AR_II=c(AR_seir2_II,AR_seir2other_II,AR_seir3_II,AR_seir3other_II)

RE_mean1=c(SEIR2_D1_RE$mean,SEIR2other_D1_RE$mean,SEIR2psi_D1_RE$mean,SEIR2psi_other_D1_RE$mean)
RE_mean2=c(SEIR2_D2_RE$mean,SEIR2other_D2_RE$mean,SEIR2psi_D2_RE$mean,SEIR2psi_other_D2_RE$mean)
RE_max1=c(SEIR2_D1_RE$max,SEIR2other_D1_RE$max,SEIR2psi_D1_RE$max,SEIR2psi_other_D1_RE$max)
RE_max2=c(SEIR2_D2_RE$max,SEIR2other_D2_RE$max,SEIR2psi_D2_RE$max,SEIR2psi_other_D2_RE$max)

# values
IC_MSE=rbind(nb_parameters,nb_obs,
             nothing,DIC,rmse_ndss,rmse_denfree,
             nothing,rmse_proj,
             R0_mean,R0_max,my_psi,
             HS,HS1,HS2,
             rep,
             nothing,AR,AR_II,
             RE_mean1,RE_mean2,RE_max1,RE_max2 )

# row names
Model=c("nb parameters","nb observations",
               "ESTIMATION SET","DIC","RMSE NDSS","RMSE DENFREE",
               "TEST SET","RMSE 2014-2015",
               "mean $R_0$","max $R_0$","$\\psi$",
               "$H_S(0)/N$ (\\%)","$H_{S1}(0)/N$ (\\%)","$H_{S2}(0)/N$ (\\%)",
               "Observation rate (\\%)",
               "Median annual incidence proportion","primary infection (\\%)","secondary infection (\\%)",
               "mean Re strain 1","mean Re strain 2",
               "max Re strain 1","max Re strain 2")
V1=c(rep("",8),rep("median (95\\%CI)",7),"",rep("median 2002-2015 (min-max)",2),rep("median (95\\%CI)",4))

IC_MSE=cbind(Model,V1,IC_MSE)
IC_MSE=IC_MSE[c("Model","V1","SEIR2_A","SEIR2_B","SEIR2psi_A","SEIR2psi_B")]
names(IC_MSE)=c("Model","","SEIR2_A","SEIR2_B","SEIR2psi_A","SEIR2psi_B")

# export table
bold <- function(x){paste0('{\\bfseries ', x, '}')}
print(xtable(IC_MSE, type = "latex"), file = file.path(dir_fig,"table_other_mode.tex"),
      include.rownames = F,hline.after = c(-1,0,6,8,11,14,15,18,22),sanitize.colnames.function = bold,
      sanitize.text.function = function(x) x,
      floating=FALSE,latex.environments=NULL,booktabs=TRUE)

