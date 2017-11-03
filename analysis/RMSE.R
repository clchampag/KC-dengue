require(data.table)
library(ggplot2)
library(gridExtra)  



setwd("/ssm")
ROOT = getwd() 
pathSEIR = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")
datapath = "data"


#RMSE on NDSS data
mse_ndss=function(path, mcmc, cut_date, burn=0 ,pattern="X_128.csv"){
  fileNames = list.files(file.path(path,mcmc), pattern = pattern, full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date","ran_incidence_ndss_obs" ,"incidence_ndss_obs","index"))
  
  fileData=list.files(file.path(path,"data"), pattern = "data.csv", full.names = T)
  data <- fread(fileData, sep = ",", select = c("date", "incidence_ndss_obs"))
  names(data)=c("date", "data")
  fileDf=fileDf[fileDf$index>=burn,]
  test=merge(fileDf, data, by="date")
  test$SE=(test$incidence_ndss_obs-test$data)^2
  test$SE_ran=(test$ran_incidence_ndss_obs-test$data)^2
  
  mse=test[, mean(SE), by = c("index")]
  mse_ran=test[, mean(SE_ran), by = c("index")]
  mse_before=test[as.Date(test$date)< cut_date, mean(SE_ran), by = c("index")]
  mse_after=test[as.Date(test$date)>= cut_date, mean(SE_ran), by = c("index")]
  
  print(c("without obs. noise",mean(sqrt(mse$V1))))
  print(c("with obs. noise",mean(sqrt(mse_ran$V1))))
  print(c("with obs. noise before ",mean(sqrt(mse_before$V1))))
  print(c("with obs. noise after",mean(sqrt(mse_after$V1))))
}


#RMSE on NDSS data (projections)
mse_ndss_proj=function(path, datapath, start_date,end_date,cut_date,mcmc){
  fileNames = list.files(file.path(path,mcmc), pattern = "X_128.csv", full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date", "ran_incidence_ndss_obs", "incidence_ndss_obs","index"))
  fileDf = subset(fileDf, as.Date(date)>"2013-12-31" & as.Date(date)<"2016-01-01")
  
  # aggregate simulations per week
  fileDf$week=as.numeric( format( as.Date(fileDf$date), "%U"))+1
  fileDf$week[fileDf$week==53]=1
  fileDf$year=as.numeric( format( as.Date(fileDf$date), "%Y"))
  fileDf$year[fileDf$year==2013]=2014
  fileDfw=fileDf[, sum(incidence_ndss_obs), by = c("year","week","index")]
  fileDfw_ran=fileDf[, sum(ran_incidence_ndss_obs), by = c("year","week","index")]
  fileDfw$week[fileDfw$week<10]=paste("0",fileDfw$week[fileDfw$week<10], sep="" )
  fileDfw$date=as.Date(paste(fileDfw$year,fileDfw$week,1,sep=""),"%Y%U%u")
  names(fileDfw)=c("year" , "week",  "index" ,"incidence_ndss_obs"  ,  "date" )
  fileDfw_ran$week[fileDfw_ran$week<10]=paste("0",fileDfw_ran$week[fileDfw_ran$week<10], sep="" )
  fileDfw_ran$date=as.Date(paste(fileDfw_ran$year,fileDfw_ran$week,1,sep=""),"%Y%U%u")
  names(fileDfw_ran)=c("year" , "week",  "index" ,"ran_incidence_ndss_obs"  ,  "date" )
  
  fileData=list.files(datapath, pattern = "data.csv", full.names = T)
  data <- fread(fileData, sep = ",", select = c("date", "incidence_ndss_obs"))
  names(data)=c("date", "data")
  data=as.data.frame(data)
  data$date=as.Date(data$date)
  fileDfw$date=as.Date(fileDfw$date)
  data=data[as.Date(data$date)<= end_date & as.Date(data$date)>= start_date,]
  test=merge(fileDfw, data, by="date")
  test_ran=merge(fileDfw_ran, data, by="date")
  
  test$SE=(test$incidence_ndss_obs-test$data)^2
  test_ran$SE_ran=(test_ran$ran_incidence_ndss_obs-test_ran$data)^2
  
  mse=test[, mean(SE), by = c("index")]
  mse_ran=test_ran[, mean(SE_ran), by = c("index")]
  mse_before=test_ran[as.Date(test_ran$date)< cut_date, mean(SE_ran), by = c("index")]
  mse_after=test_ran[as.Date(test_ran$date)>= cut_date, mean(SE_ran), by = c("index")]

  print(c("without obs. noise",mean(sqrt(mse$V1))))
  print(c("with obs. noise",mean(sqrt(mse_ran$V1))))
  print(c("with obs. noise before ",mean(sqrt(mse_before$V1))))
  print(c("with obs. noise after",mean(sqrt(mse_after$V1))))
}


SEIR=mse_ndss(pathSEIR,mcmc = "mcmc",cut_date = "2003-01-01")
Laneri=mse_ndss(pathLaneri,mcmc = "mcmc", cut_date = "2003-01-01")
Pandey=mse_ndss(pathPandey,mcmc = "mcmc", cut_date = "2003-01-01")
SEIAR=mse_ndss(pathSEIAR,mcmc = "mcmc", cut_date = "2003-01-01")
SEIR2=mse_ndss(pathSEIR2,mcmc = "mcmc", cut_date = "2003-01-01")
SEIR2_3=mse_ndss(pathSEIR2_psi,mcmc = "mcmc", cut_date = "2003-01-01")


SEIR=mse_ndss_proj(pathSEIR,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
Laneri=mse_ndss_proj(pathLaneri,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
Pandey=mse_ndss_proj(pathPandey,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIAR=mse_ndss_proj(pathSEIAR,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIR2=mse_ndss_proj(pathSEIR2,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")
SEIR2_3=mse_ndss_proj(pathSEIR2_psi,datapath , cut_date = "2014-12-31",end_date = "2015-12-31",start_date = "2014-01-01", mcmc="mcmc")

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
}

mse_denfree(pathSEIR,mcmc="mcmc")
mse_denfree(pathLaneri,mcmc="mcmc")
mse_denfree(pathPandey,mcmc="mcmc")
mse_denfree(pathSEIAR,mcmc="mcmc")
mse_SEIR2(pathSEIR2,mcmc="mcmc")
mse_SEIR2(pathSEIR2_psi,mcmc="mcmc")
