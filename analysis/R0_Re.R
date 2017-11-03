require(data.table)
setwd("/ssm")

ROOT = getwd() 
pathSEIR = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")
N=161391

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
  medianR=apply(var_R0, 2, median)
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  
  print(c("median",quantile(medianR   , probs=c(0.025,0.5,0.975))))
  print(c("mean",quantile(meanR   , probs=c(0.025,0.5,0.975))))
  print(c("max",quantile(maxR   , probs=c(0.025,0.5,0.975))))
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
  medianR=apply(var_R0, 2, median)
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  minR=apply(var_R0, 2, min)
  
  print(c("median",quantile(medianR   , probs=c(0.025,0.5,0.975))))
  print(c("mean",quantile(meanR   , probs=c(0.025,0.5,0.975))))
  print(c("max",quantile(maxR   , probs=c(0.025,0.5,0.975))))
  print(c("min",quantile(minR   , probs=c(0.025,0.5,0.975))))
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
  medianR=apply(var_R0, 2, median)
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  minR=apply(var_R0, 2, min)
  
  print(c("median",quantile(medianR   , probs=c(0.025,0.5,0.975))))
  print(c("mean",quantile(meanR   , probs=c(0.025,0.5,0.975))))
  print(c("max",quantile(maxR   , probs=c(0.025,0.5,0.975))))
  print(c("min",quantile(minR   , probs=c(0.025,0.5,0.975))))
  
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
  medianR=apply(var_R0, 2, median)
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  minR=apply(var_R0, 2, min)
  
  print(c("median",quantile(medianR   , probs=c(0.025,0.5,0.975))))
  print(c("mean",quantile(meanR   , probs=c(0.025,0.5,0.975))))
  print(c("max",quantile(maxR   , probs=c(0.025,0.5,0.975))))
  print(c("min",quantile(minR   , probs=c(0.025,0.5,0.975))))
}  

# Effective reproduction umber (two strain models with interaction (psi))
Reff_strains_psi = function(path,N,S,S1, I2, I12,name){
  fileNames = list.files(file.path(path,"mcmc"), pattern = paste0("X_128.csv"), full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date",S,S1,I2,I12,"index"))
  names(fileDf)=c("date","S","S1","I2","I12","index")
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
    var_R0[i,] = ((test$I2[test$x==xx]/(test$I2[test$x==xx]+test$I12[test$x==xx]))+(test$psi[test$x==xx]*test$I12[test$x==xx]/(test$I2[test$x==xx]+test$I12[test$x==xx])))*((test$S[test$x==xx]+test$S1[test$x==xx])/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  medianR=apply(var_R0, 2, median)
  meanR=apply(var_R0, 2, mean)
  maxR=apply(var_R0, 2, max)
  minR=apply(var_R0, 2, min)
  
  print(c("median",quantile(medianR   , probs=c(0.025,0.5,0.975))))
  print(c("mean",quantile(meanR   , probs=c(0.025,0.5,0.975))))
  print(c("max",quantile(maxR   , probs=c(0.025,0.5,0.975))))
  print(c("min",quantile(minR   , probs=c(0.025,0.5,0.975))))
}  

R0_S0(pathSEIR)
R0_S0(pathLaneri)
R0_S0(pathPandey)
R0_S0(pathSEIAR)
R0_S0(pathSEIR2)
R0_S0(pathSEIR2_psi)


reff_onestrain(pathSEIR,N,"S", name="seir")
reff_onestrain(pathLaneri,N,"S",name="laneri")
ReffPandey=Reff_pandey(pathPandey,N,name="pandey")
ReffSEIAR=reff_onestrain(pathSEIAR,N,"S",name="seiar")

Reff_strains(pathSEIR2,N,"S","S2",name="seir2_D1")
Reff_strains(pathSEIR2,N,"S","S1",name="seir2_D2")

Reff_strains_psi(pathSEIR2_psi,N,"S","S2","I1","I21",name="seir2psi_D1")
Reff_strains_psi(pathSEIR2_psi,N,"S","S1","I2","I12",name="seir2psi_D2")
