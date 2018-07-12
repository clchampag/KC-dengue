library(gridExtra)  
require(data.table)
library(ggplot2)


dir_fig = "."
setwd("/ssm")


#setwd("/ssm")
ROOT = getwd() 
pathSEIR = paste(ROOT,"/SEIR", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")
N=161391

proj=function(path, end_date, state,to_sum,mcmc,X){
  fileNames = list.files(file.path(path,mcmc), pattern = X, full.names = T)
  fileDf = fread(fileNames, sep = ",", select = c("date", state,"index"))
  names(fileDf)=c("date",state,"index")
  
  median=fileDf[, median(eval(parse(text = to_sum))), by = c("date")]
  q975=fileDf[, quantile(eval(parse(text = to_sum)), probs=c(0.975)), by = c("date")]
  q025=fileDf[, quantile(eval(parse(text = to_sum)), probs=c(0.025)), by = c("date")]
  
  tot=data.frame(date=as.Date(median$date), median=median$V1, q975=q975$V1, q025=q025$V1)
  tot=tot[as.Date(tot$date)<end_date,]
  
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  ggplot(tot,aes_string(x='date',y ="median")) +
    geom_ribbon(aes_string(ymin= "q025", ymax="q975"),fill="skyblue",alpha=0.55) +
    geom_line(color="dodgerblue4",size=0.7)+ theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))+
    theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
    theme(axis.title = element_text(size=16)) +ylab("") + xlab("") +
    theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))
  
}

# Plot each model separately
SEIR2_I1=proj(pathSEIR2_psi,end_date = "2013-12-31",state=c("I1","I21"),to_sum="I1+I21",mcmc="mcmc",X="X_128.csv")
SEIR2_I2=proj(pathSEIR2_psi,end_date = "2013-12-31",state=c("I2","I12"),to_sum="I2+I12",mcmc="mcmc",X="X_128.csv")
SEIR2_S1=proj(pathSEIR2_psi,end_date = "2013-12-31",state=c("S1"),to_sum="S1/N",mcmc="mcmc",X="X_128.csv")
SEIR2_S2=proj(pathSEIR2_psi,end_date = "2013-12-31",state=c("S2"),to_sum="S2/N",mcmc="mcmc",X="X_128.csv")
SEIR2_S=proj(pathSEIR2_psi,end_date = "2013-12-31",state=c("S"),to_sum="S/N",mcmc="mcmc",X="X_128.csv")
SEIR_S=proj(pathSEIR,end_date = "2013-12-31",state=c("S"),to_sum="S/N",mcmc="mcmc",X="X_128.csv")

#Graphical settings
SEIR2_I1b=SEIR2_I1+ ylim(0,1000)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("A") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_I2b=SEIR2_I2+ ylim(0,1000)+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("B") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR_Sb=SEIR_S+ ylim(0,0.5)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("C") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S1b=SEIR2_S1+ ylim(0,0.5)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("D") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S2b=SEIR2_S2+ ylim(0,0.5)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("E") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_Sb=SEIR2_S+ ylim(0,0.5)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("F") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))

#Combine plots
hlay <- rbind(c(1,2,3),
              c(4,5,6))

plot=grid.arrange(grobs=list(SEIR2_I1b,SEIR2_I2b,SEIR_Sb,SEIR2_S1b,SEIR2_S2b,SEIR2_Sb), layout_matrix=hlay)
ggsave(file.path(dir_fig,"susceptibles.pdf"),plot, width=40, height=30)


