require(data.table)
library(dplyr)
library(ggplot2)


dir_fig = "."
setwd("/ssm")
ROOT = getwd() 
pathSEIR = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")
N=161391



epidemio_attack_rate=function(path, S, to_sum_S, incidence, to_sum,mcmc){
  fileNames = list.files(file.path(path,mcmc), pattern = "X_128.csv*", full.names = T)
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
  q25_attack_rate=tot[, quantile(attack_rate, probs=c(0.025)), by = c("year")]
  names(q25_attack_rate)=c("date","q025")
  q95_attack_rate=tot[, quantile(attack_rate, probs=c(0.975)), by = c("year")]
  names(q95_attack_rate)=c("date","q975")
  
  AR=merge(med_attack_rate, q25_attack_rate, by=c("date"))
  AR=merge(AR, q95_attack_rate, by=c("date"))
  
  return(AR)
}


###### ATTACK RATE
AR_seir=epidemio_attack_rate(pathSEIR,S="S",to_sum_S = "S",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_laneri=epidemio_attack_rate(pathLaneri,S="S",to_sum_S = "S",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_pandey=epidemio_attack_rate(pathPandey,S="Hs",to_sum_S = "Hs",incidence = "incidence", to_sum="incidence",mcmc="mcmc")
AR_seiar=epidemio_attack_rate(pathSEIAR,S="S",to_sum_S = "S",incidence=c("incidence_I", "incidence_H","incidence_asympto"),to_sum = "incidence_I+incidence_H+incidence_asympto",mcmc="mcmc")
AR_seir2_1=epidemio_attack_rate(pathSEIR2,S=c("S","S2"),to_sum_S="S+S2",incidence=c("inc_I1","inc_I21"), to_sum="inc_I1+inc_I21",mcmc="mcmc")
AR_seir2_2=epidemio_attack_rate(pathSEIR2,S=c("S","S1"),to_sum_S="S+S1",incidence=c("inc_I2","inc_I12"), to_sum="inc_I2+inc_I12",mcmc="mcmc")
AR_seir3_1=epidemio_attack_rate(pathSEIR2_psi,S=c("S","S2"),to_sum_S="S+S2",incidence=c("inc_I1","inc_I21"), to_sum="inc_I1+inc_I21",mcmc="mcmc")
AR_seir3_2=epidemio_attack_rate(pathSEIR2_psi,S=c("S","S1"),to_sum_S="S+S1",incidence=c("inc_I2","inc_I12"), to_sum="inc_I2+inc_I12",mcmc="mcmc")
AR_seir2_I=epidemio_attack_rate(pathSEIR2,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2",mcmc="mcmc")
AR_seir3_I=epidemio_attack_rate(pathSEIR2_psi,S=c("S"),to_sum_S="S",incidence=c("inc_I1","inc_I2"), to_sum="inc_I1+inc_I2",mcmc="mcmc")
AR_seir2_II=epidemio_attack_rate(pathSEIR2,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I12","inc_I21"), to_sum="inc_I12+inc_I21",mcmc="mcmc")
AR_seir3_II=epidemio_attack_rate(pathSEIR2_psi,S=c("S1","S2"),to_sum_S="S1+S2",incidence=c("inc_I21","inc_I12"), to_sum="inc_I12+inc_I21",mcmc="mcmc")

summary(AR_seir$med)
summary(AR_laneri$med)
summary(AR_pandey$med)
summary(AR_seiar$med)
summary(AR_seir2_1$med)
summary(AR_seir2_2$med)
summary(AR_seir3_1$med)
summary(AR_seir3_2$med)
summary(AR_seir2_I$med)
summary(AR_seir3_I$med)
summary(AR_seir2_II$med)
summary(AR_seir3_II$med)

names(AR_seir)=c("date", "med_seir",  "q025_seir" ,"q975_seir")
names(AR_laneri)=c("date", "med_lan",  "q025_lan" ,"q975_lan" )
names(AR_pandey)=c("date", "med_pandey",  "q025_pandey" ,"q975_pandey")
names(AR_seiar)=c("date", "med_seiar",  "q025_seiar" ,"q975_seiar" )
names(AR_seir2_I)=c("date", "med_seir2I",  "q025_seir2I" ,"q975_seir2I" )
names(AR_seir3_I)=c("date", "med_seir3I",  "q025_seir3I" ,"q975_seir3I" )
names(AR_seir2_II)=c("date", "med_seir2II",  "q025_seir2II" ,"q975_seir2II" )
names(AR_seir3_II)=c("date", "med_seir3II",  "q025_seir3II" ,"q975_seir3II" )
AR_all=merge(AR_seir, AR_laneri)
AR_all=merge(AR_all, AR_pandey)
AR_all=merge(AR_all, AR_seiar)
AR_all=merge(AR_all, AR_seir2_I)
AR_all=merge(AR_all, AR_seir3_I)
AR_all=merge(AR_all, AR_seir2_II)
AR_all=merge(AR_all, AR_seir3_II)

myvars=c("date","med_seir","med_lan","med_pandey","med_seiar","med_seir2I", "med_seir3I")
AR_all=as.data.frame(AR_all)
mdata <- melt(AR_all[myvars], id=c("date")) 

ar=ggplot(mdata,aes_string(x='date')) +geom_point(aes_string( y="value", group="variable", color="variable"),size=4)+
  scale_colour_manual(values=c("skyblue", "dodgerblue", "darkblue","forestgreen","red","orange"),
                      name="",
                      labels=c("SEIR", "Laneri", "Pandey","SEIAR","SEIR2","SEIR2psi"))+
  geom_errorbar(data=AR_seir,aes_string(ymin= "q025_seir", ymax="q975_seir"),color="skyblue",alpha=1, size=1.5)+
  geom_errorbar(data=AR_laneri,aes_string(ymin= "q025_lan", ymax="q975_lan"),color="dodgerblue",alpha=1, size=1.5) +
  geom_errorbar(data=AR_pandey,aes_string(ymin= "q025_pandey", ymax="q975_pandey"),color="darkblue",alpha=1, size=1.5) +
  geom_errorbar(data=AR_seiar,aes_string(ymin= "q025_seiar", ymax="q975_seiar"),color="forestgreen",alpha=1, size=1.5) +
  
  geom_errorbar(data=AR_seir2_I,aes_string(ymin= "q025_seir2I", ymax="q975_seir2I"),color="red",alpha=1, size=1.5) +
  geom_errorbar(data=AR_seir3_I,aes_string(ymin= "q025_seir3I", ymax="q975_seir3I"),color="orange",alpha=1, size=1.5) +
  
  theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
  theme(axis.title = element_text(size=16)) +ylab("") + xlab("") +
  theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))+
  theme(legend.key.size = unit(2, "cm"), legend.text = element_text(size=30))

ggsave(file.path(dir_fig,"inc.pdf"),ar, width=18, height=6)
