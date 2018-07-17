require(data.table)
library(dplyr)


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
N=161391

ci_trace=function(path,parameter, to_sum, my_digits=0){
  fileNames = list.files(file.path(path,"mcmc"), pattern ="trace", full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c(parameter))
  fileDf=mutate(fileDf,param=eval(parse(text = to_sum)))
  
  my_ci=quantile(fileDf$param, probs = c(0.5,0.025,0.975))
  my_ci=round(my_ci,digits = my_digits)
  return(paste0(my_ci[1]," (",my_ci[2],"-",my_ci[3],")"))
}  


# Basic reproduction number
R0_S0 = function(path){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("trace_1.csv"), full.names = T)
  fileTrace = fread(fileNames[1], sep = ",")
  date=seq.Date(as.Date("2002-01-14"),as.Date("2003-01-14"),by="day")
  var_R0=matrix(nrow = length(date), ncol=dim(fileTrace)[1]+1)
  var_R0[,1]=seq(1:length(date))
  for(i in (1:length(date))){
    var_R0[i,-1] = fileTrace$R0*(1 + fileTrace$beta1*sin(2*pi*(var_R0[i,1]/365+fileTrace$phase)))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_ci_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_ci_mean=round(my_ci_mean,digits = 2)
  my_ci_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_ci_max=round(my_ci_max,digits = 2)
  return(list(mean=paste0(my_ci_mean[1]," (",my_ci_mean[2],"-",my_ci_mean[3],")"),
              max=paste0(my_ci_max[1]," (",my_ci_max[2],"-",my_ci_max[3],")")))
}  

# Effective reproduction umber (one strain models)
reff_onestrain = function(path,N,S,name){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date", S,"index"))
  names(fileDf)=c("date","S","index")
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("trace_1.csv"), full.names = T)
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
    var_R0[i,] = (test$S[test$x==xx]/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_ci_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_ci_mean=round(my_ci_mean,digits = 2)
  my_ci_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_ci_max=round(my_ci_max,digits = 2)
  return(list(mean=paste0(my_ci_mean[1]," (",my_ci_mean[2],"-",my_ci_mean[3],")"),
              max=paste0(my_ci_max[1]," (",my_ci_max[2],"-",my_ci_max[3],")")))
}  

# Effective reproduction umber (Pandey model)
Reff_pandey = function(path,N,name){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date","Hs","Vs","index"))
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("trace_1.csv"), full.names = T)
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
    var_R0[i,] = (test$Hs[test$x==xx]/N)*(test$Vs[test$x==xx])*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_ci_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_ci_mean=round(my_ci_mean,digits = 2)
  my_ci_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_ci_max=round(my_ci_max,digits = 2)
  return(list(mean=paste0(my_ci_mean[1]," (",my_ci_mean[2],"-",my_ci_mean[3],")"),
              max=paste0(my_ci_max[1]," (",my_ci_max[2],"-",my_ci_max[3],")")))
}  

# Effective reproduction number (two strain model)
Reff_strains = function(path,N,S,S1,name){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",S,S1,"index"))
  names(fileDf)=c("date","S","S1","index")
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("trace_1.csv"), full.names = T)
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
  
  my_ci_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_ci_mean=round(my_ci_mean,digits = 2)
  my_ci_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_ci_max=round(my_ci_max,digits = 2)
  return(list(mean=paste0(my_ci_mean[1]," (",my_ci_mean[2],"-",my_ci_mean[3],")"),
              max=paste0(my_ci_max[1]," (",my_ci_max[2],"-",my_ci_max[3],")")))
}  

# Effective reproduction umber (two strain models with interaction (psi))
Reff_strains_psi = function(path,N,S,S1, name){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",S,S1,"index"))
  names(fileDf)=c("date","S","S1","index")
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("trace_1.csv"), full.names = T)
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
    var_R0[i,] = ((test$S[test$x==xx]+test$psi[test$x==xx]*test$S1[test$x==xx])/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  my_ci_mean=quantile(meanR, probs=c(0.5,0.025,0.975))
  my_ci_mean=round(my_ci_mean,digits = 2)
  my_ci_max=quantile(maxR, probs=c(0.5,0.025,0.975))
  my_ci_max=round(my_ci_max,digits = 2)
  return(list(mean=paste0(my_ci_mean[1]," (",my_ci_mean[2],"-",my_ci_mean[3],")"),
              max=paste0(my_ci_max[1]," (",my_ci_max[2],"-",my_ci_max[3],")")))
}  


epidemio_attack_rate=function(path, S, to_sum_S, incidence, to_sum,mcmc){
  fileNames = list.files(file.path(path,mcmc), pattern = "X_128.csv", full.names = T)
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
SEIR_R0=R0_S0(pathSEIR)
Laneri_R0=R0_S0(pathLaneri)
Pandey_R0=R0_S0(pathPandey)
SEIAR_R0=R0_S0(pathSEIAR)
SEIR2_R0=R0_S0(pathSEIR2)
SEIR2psi_R0=R0_S0(pathSEIR2_psi)

# RE
SEIR_RE=reff_onestrain(pathSEIR,N,"S", name="seir")
Laneri_RE=reff_onestrain(pathLaneri,N,"S",name="laneri")
Pandey_RE=ReffPandey=Reff_pandey(pathPandey,N,name="pandey")
SEIAR_RE=ReffSEIAR=reff_onestrain(pathSEIAR,N,"S",name="seiar")
SEIR2_D1_RE=Reff_strains(pathSEIR2,N,"S","S2",name="seir2_D1")
SEIR2_D2_RE=Reff_strains(pathSEIR2,N,"S","S1",name="seir2_D2")
SEIR2psi_D1_RE=Reff_strains_psi(path=pathSEIR2_psi,N,S="S",S1="S2",name="seir2psi_D1")
SEIR2psi_D2_RE=Reff_strains_psi(pathSEIR2_psi,N,"S","S1",name="seir2psi_D2")



# Susceptibles
SEIR_S=ci_trace(pathSEIR,"pr_S",to_sum = "pr_S*100")
Laneri_S=ci_trace(pathLaneri,"pr_S",to_sum = "pr_S*100")
Pandey_S=ci_trace(pathPandey,"pr_Hs",to_sum = "pr_Hs*100")
SEIAR_S=ci_trace(pathSEIAR,"pr_S",to_sum = "pr_S*100")
SEIR2_S=ci_trace(pathSEIR2,"pr_S",to_sum = "pr_S*100")
SEIR2psi_S=ci_trace(pathSEIR2_psi,"pr_S",to_sum = "pr_S*100")

SEIR2_S1=ci_trace(pathSEIR2,"pr_S1",to_sum = "pr_S1*100")
SEIR2psi_S1=ci_trace(pathSEIR2_psi,"pr_S1",to_sum = "pr_S1*100")
SEIR2_S2=ci_trace(pathSEIR2,"pr_S2",to_sum = "pr_S2*100")
SEIR2psi_S2=ci_trace(pathSEIR2_psi,"pr_S2",to_sum = "pr_S2*100")


SEIR_r=ci_trace(pathSEIR,"rep_ndss","rep_ndss*100")
Laneri_r=ci_trace(pathLaneri,"rep_ndss","rep_ndss*100")
Pandey_r=ci_trace(pathPandey,"rep_ndss","rep_ndss*100")
SEIAR_rh=ci_trace(pathSEIAR,"rho_h","rho_h*100")
SEIAR_ra=ci_trace(pathSEIAR,"rho_a","rho_a*100")
SEIR2_r=ci_trace(pathSEIR2,"rep_ndss","rep_ndss*100")
SEIR2psi_r=ci_trace(pathSEIR2_psi,"rep_ndss","rep_ndss*100")

psi=ci_trace(pathSEIR2_psi,"psi","psi", my_digits = 2)

SEIR_disp=ci_trace(pathSEIR,"disp","disp", my_digits = 2)
Laneri_disp=ci_trace(pathLaneri,"disp","disp", my_digits = 2)
Pandey_disp=ci_trace(pathPandey,"disp","disp", my_digits = 2)
SEIAR_disp=ci_trace(pathSEIAR,"disp","disp", my_digits = 2)
SEIR2_disp=ci_trace(pathSEIR2,"disp","disp", my_digits = 2)
SEIR2psi_disp=ci_trace(pathSEIR2_psi,"disp","disp", my_digits = 2)



AR_seir=epidemio_attack_rate(pathSEIR,S="S",to_sum_S = "S",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_laneri=epidemio_attack_rate(pathLaneri,S="S",to_sum_S = "S",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_pandey=epidemio_attack_rate(pathPandey,S="Hs",to_sum_S = "Hs",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_seiar=epidemio_attack_rate(pathSEIAR,S="S",to_sum_S = "S",incidence=c("incidence_I", "incidence_H","incidence_asympto"),to_sum = "incidence_I+incidence_H+incidence_asympto",mcmc="mcmc")
AR_seir2_I=epidemio_attack_rate(pathSEIR2,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2",mcmc="mcmc")
AR_seir3_I=epidemio_attack_rate(pathSEIR2_psi,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2",mcmc="mcmc")
AR_seir2_II=epidemio_attack_rate(pathSEIR2,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I12","inc_I21"), to_sum="inc_I12+inc_I21",mcmc="mcmc")
AR_seir3_II=epidemio_attack_rate(pathSEIR2_psi,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I21","inc_I12"), to_sum="inc_I12+inc_I21",mcmc="mcmc")


# create latex table
library(xtable)
nothing=c("","","","","","")
R0_mean=c(SEIR_R0$mean,Laneri_R0$mean, Pandey_R0$mean, SEIAR_R0$mean,SEIR2_R0$mean, SEIR2psi_R0$mean)
names(R0_mean)=c("SEIR","Laneri","Pandey","SEIAR","SEIR2","SEIR2 psi")
R0_max=c(SEIR_R0$max,Laneri_R0$max, Pandey_R0$max, SEIAR_R0$max,SEIR2_R0$max, SEIR2psi_R0$max)
my_psi=c("","","","","",psi)

HS=c(SEIR_S,Laneri_S, Pandey_S, SEIAR_S,SEIR2_S, SEIR2psi_S)
HS1=c("","","","",SEIR2_S1, SEIR2psi_S1)
HS2=c("","","","",SEIR2_S2, SEIR2psi_S2)

rep=c(SEIR_r,Laneri_r, Pandey_r,"",SEIR2_r, SEIR2psi_r)
rho_h=c("","","",SEIAR_rh,"","")
rho_a=c("","","",SEIAR_ra,"","")

disp=c(SEIR_disp,Laneri_disp, Pandey_disp,SEIAR_disp,SEIR_disp, SEIR2psi_disp)

AR=c(AR_seir,AR_laneri,AR_pandey,AR_seiar,AR_seir2_I,AR_seir3_I)
AR_II=c("","","","",AR_seir2_II,AR_seir3_II)

RE_mean=c(SEIR_RE$mean,Laneri_RE$mean, Pandey_RE$mean, SEIAR_RE$mean,"","")
RE_mean1=c("","","","",SEIR2_D1_RE$mean,SEIR2psi_D1_RE$mean)
RE_mean2=c("","","","",SEIR2_D2_RE$mean,SEIR2psi_D2_RE$mean)
RE_max=c(SEIR_RE$max,Laneri_RE$max, Pandey_RE$max, SEIAR_RE$max,"","")
RE_max1=c("","","","",SEIR2_D1_RE$max,SEIR2psi_D1_RE$max)
RE_max2=c("","","","",SEIR2_D2_RE$max,SEIR2psi_D2_RE$max)

# values
EPIDEMIO=rbind(R0_mean,R0_max,my_psi,
               HS,HS1,HS2,
               rep,rho_h,rho_a,disp,
               nothing,AR,AR_II,
               RE_mean,RE_mean1,RE_mean2,RE_max,RE_max1,RE_max2 )

# row names
Model=c("mean $R_0$","max $R_0$","$\\psi$",
        "$H_S(0)/N$ (\\%)","$H_{S1}(0)/N$ (\\%)","$H_{S2}(0)/N$ (\\%)",
        "Observation rate (\\%)","Hospitalized (\\%)","Asymptomatic (\\%)","Over-dispersion",
        "Median annual incidence proportion","primary infection (\\%)","secondary infection (\\%)",
        "mean Re","mean Re strain 1","mean Re strain 2",
        "max Re","max Re strain 1","max Re strain 2")
V1=c(rep("median (95\\%CI)",10),"",rep("median 2002-2015 (min-max)",2),rep("median (95\\%CI)",6))

EPIDEMIO=cbind(Model,V1,EPIDEMIO)
names(EPIDEMIO)=c("Model","","SEIR","Laneri","Pandey","SEIAR","SEIR2","SEIR2psi")
# export table
bold <- function(x){paste0('{\\bfseries ', x, '}')}
print(xtable(EPIDEMIO, type = "latex"), file = file.path(dir_fig,"epidemio.tex"),
      include.rownames = F,hline.after = c(-1,0,3,6,10,13,19),sanitize.colnames.function = bold,
      sanitize.text.function = function(x) x,
      floating=FALSE,latex.environments=NULL,booktabs=TRUE)

