require(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)  

dir_fig = "./supp"
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
  var_R0=as.data.frame(var_R0)
  quantiles=as.data.frame(t(apply(var_R0, FUN = quantile, MARGIN = 1, probs=c(0.025,0.5,0.975))))
  names(quantiles)=c(paste("q025",name,sep = "_"),paste("q50",name,sep = "_"),paste("q975",name,sep = "_"))
  quantiles$date=test$date[test$index==0]
  return(quantiles)
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
  var_R0=as.data.frame(var_R0)
  quantiles=as.data.frame(t(apply(var_R0, FUN = quantile, MARGIN = 1, probs=c(0.025,0.5,0.975))))
  names(quantiles)=c(paste("q025",name,sep = "_"),paste("q50",name,sep = "_"),paste("q975",name,sep = "_"))
  quantiles$date=test$date[test$index==0]
  return(quantiles)
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
  var_R0=as.data.frame(var_R0)
  quantiles=as.data.frame(t(apply(var_R0, FUN = quantile, MARGIN = 1, probs=c(0.025,0.5,0.975))))
  names(quantiles)=c(paste("q025",name,sep = "_"),paste("q50",name,sep = "_"),paste("q975",name,sep = "_"))
  quantiles$date=test$date[test$index==0]
  return(quantiles)
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
    var_R0[i,] = ((test$S[test$x==xx]+test$psi[test$x==xx]*test$S1[test$x==xx])/N)*test$R0[test$x==xx]*(1 + test$beta1[test$x==xx]*sin(2*pi*(xx/365+test$phase[test$x==xx])))
  }
  var_R0=as.data.frame(var_R0)
  quantiles=as.data.frame(t(apply(var_R0, FUN = quantile, MARGIN = 1, probs=c(0.025,0.5,0.975))))
  names(quantiles)=c(paste("q025",name,sep = "_"),paste("q50",name,sep = "_"),paste("q975",name,sep = "_"))
  quantiles$date=test$date[test$index==0]
  return(quantiles)
}  


# RE
SEIR_RE=reff_onestrain(pathSEIR,N,"S", name="seir")
Laneri_RE=reff_onestrain(pathLaneri,N,"S",name="laneri")
Pandey_RE=ReffPandey=Reff_pandey(pathPandey,N,name="pandey")
SEIAR_RE=ReffSEIAR=reff_onestrain(pathSEIAR,N,"S",name="seiar")
SEIR2_D1_RE=Reff_strains(pathSEIR2,N,"S","S2",name="seir2_D1")
SEIR2_D2_RE=Reff_strains(pathSEIR2,N,"S","S1",name="seir2_D2")
SEIR2psi_D1_RE=Reff_strains_psi(path=pathSEIR2_psi,N,S="S",S1="S2",I2="I1",I12="I21",name="seir2psi_D1")
SEIR2psi_D2_RE=Reff_strains_psi(pathSEIR2_psi,N,"S","S1","I2","I12",name="seir2psi_D2")

RE_all=merge(SEIR_RE, Laneri_RE)
RE_all=merge(RE_all, SEIAR_RE)
RE_all=merge(RE_all, Pandey_RE)
RE_all=merge(RE_all, SEIR2_D1_RE)
RE_all=merge(RE_all, SEIR2_D2_RE)
RE_all=merge(RE_all, SEIR2psi_D1_RE)
RE_all=merge(RE_all, SEIR2psi_D2_RE)

# ONE STRAIN MODELS
myvars=c("date","q50_seir","q50_laneri","q50_pandey","q50_seiar")
mdata <- melt(RE_all[myvars], id=c("date")) 
mdata$date=as.Date(mdata$date)

re=ggplot(mdata,aes_string(x='date')) +geom_line(aes_string( y="value", group="variable", color="variable"),size=3)+
  scale_colour_manual(values=c("skyblue", "forestgreen", "darkblue","dodgerblue"),
                      name="",
                      labels=c("SEIR", "Laneri", "Pandey","SEIAR"))+
  
  theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
  theme(axis.title = element_text(size=16)) +ylab("") + xlab("") + ylim(0,2.5)+
  ggtitle('A') +theme(plot.title=element_text( face="bold", size=45,hjust =0))+
  theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))+
  theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size=30))

# SEIR2 model

myvars2=c("date","q50_seir2_D1","q50_seir2_D2")
mdata2 <- melt(RE_all[myvars2], id=c("date")) 
mdata2$date=as.Date(mdata2$date)



re2=ggplot(mdata2,aes_string(x='date')) +geom_line(aes_string( y="value", group="variable", color="variable"),size=3)+
  scale_colour_manual(values=c("orange","dodgerblue"),
                      name="",
                      labels=c("Strain 1","Strain 2"))+
  
  theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
  theme(axis.title = element_text(size=16)) +ylab("") + xlab("") + ylim(0,2.5)+
  ggtitle('B') +theme(plot.title=element_text( face="bold", size=45,hjust =0))+
  theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))+
  theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size=30))

myvars2=c("date","q50_seir2_D1", "q50_seir2psi_D1","q50_seir2_D2", "q50_seir2psi_D2")
mdata2 <- melt(RE_all[myvars2], id=c("date")) 
mdata2$date=as.Date(mdata2$date)

myvars3=c("date","q50_seir2psi_D1","q50_seir2psi_D2")
mdata3 <- melt(RE_all[myvars3], id=c("date")) 
mdata3$date=as.Date(mdata3$date)

# SEIR2psi model

re3=ggplot(mdata3,aes_string(x='date')) +geom_line(aes_string( y="value", group="variable", color="variable"),size=3)+
  scale_colour_manual(values=c("orange","dodgerblue"),
                      name="",
                      labels=c("Strain 1","Strain 2"))+
  
  theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
  theme(axis.title = element_text(size=16)) +ylab("") + xlab("") + ylim(0,2.5)+
  ggtitle('C') +theme(plot.title=element_text( face="bold", size=45,hjust =0))+
  theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))+
  theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size=30))


hlay <- rbind(c(1),
              c(2),
              c(3))

plot_re=grid.arrange(grobs=list(re,re2,re3),  layout_matrix=hlay)


ggsave(file.path(dir_fig,"re.pdf"),plot_re, width=18, height=18)

