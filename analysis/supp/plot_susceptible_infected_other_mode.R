library(gridExtra)  
require(data.table)
library(ggplot2)

dir_fig = "./supp"
setwd("/ssm")


#setwd("/ssm")
ROOT = getwd() 
pathSEIR2 = paste(ROOT,"/SEIR2/mcmc", sep="")
pathSEIR2_paramB = paste(ROOT,"/SEIR2/mcmc-paramB", sep="")

pathSEIR2_psi = paste(ROOT,"/SEIR2-psi/mcmc", sep="")
pathSEIR2_psi_paramB = paste(ROOT,"/SEIR2-psi/mcmc-paramB", sep="")

N=161391


proj=function(path1,path2, end_date, state,to_sum,X){
  fileName1 = list.files(path1, pattern = X, full.names = T)
  fileDf1 = fread(fileName1, sep = ",", select = c("date", state,"index"))
  names(fileDf1)=c("date",state,"index")
  fileName2 = list.files(path2, pattern = X, full.names = T)
  fileDf2 = fread(fileName2, sep = ",", select = c("date", state,"index"))
  names(fileDf2)=c("date",state,"index")

  median1=fileDf1[, median(eval(parse(text = to_sum))), by = c("date")]
  q975_1=fileDf1[, quantile(eval(parse(text = to_sum)), probs=c(0.975)), by = c("date")]
  q025_1=fileDf1[, quantile(eval(parse(text = to_sum)), probs=c(0.025)), by = c("date")]
  median2=fileDf2[, median(eval(parse(text = to_sum))), by = c("date")]
  q975_2=fileDf2[, quantile(eval(parse(text = to_sum)), probs=c(0.975)), by = c("date")]
  q025_2=fileDf2[, quantile(eval(parse(text = to_sum)), probs=c(0.025)), by = c("date")]

  tot=data.frame(date=as.Date(median1$date), median1=median1$V1, q975_1=q975_1$V1, q025_1=q025_1$V1,
                 median2=median2$V1, q975_2=q975_2$V1, q025_2=q025_2$V1)
  tot=tot[as.Date(tot$date)<end_date,]
  
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  ggplot(tot,aes_string(x='date',y ="median1")) +
    geom_ribbon(aes_string(ymin= "q025_1", ymax="q975_1"),fill="skyblue",alpha=0.55) +
    geom_line(color="dodgerblue4",size=0.7) +
    geom_ribbon(aes_string(ymin= "q025_2", ymax="q975_2"),fill="red",alpha=0.55) +
    geom_line(aes_string(y ="median2"),color="red",size=0.7)+ theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))+
    theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
    theme(axis.title = element_text(size=16)) +ylab("") + xlab("") +
    theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))
  
}



# Plot each model separately
SEIR2_I1=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("I1","I21"),to_sum="I1+I21",X="X_128.csv")
SEIR2_I2=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("I2","I12"),to_sum="I2+I12",X="X_128.csv")
SEIR2_I=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("I2","I12","I1","I21"),to_sum="I2+I12+I1+I21",X="X_128.csv")
SEIR2_S1=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("S1"),to_sum="S1/N",X="X_128.csv")
SEIR2_S2=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("S2"),to_sum="S2/N",X="X_128.csv")
SEIR2_S=proj(path1=pathSEIR2, path2=pathSEIR2_paramB,end_date = "2013-12-31",state=c("S"),to_sum="S/N",X="X_128.csv")




#Graphical settings
SEIR2_I1b=SEIR2_I1+ ylim(0,1000)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("A") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_I2b=SEIR2_I2+ ylim(0,1000)+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("B") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_Ib=SEIR2_I+ ylim(0,1000)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("C") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S1b=SEIR2_S1+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("D") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S2b=SEIR2_S2+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("E") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_Sb=SEIR2_S+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("F") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))

#Combine plots
hlay <- rbind(c(1,2,3),
              c(4,5,6))

plot=grid.arrange(grobs=list(SEIR2_I1b,SEIR2_I2b,SEIR2_Ib,SEIR2_S1b,SEIR2_S2b,SEIR2_Sb), layout_matrix=hlay)
ggsave(file.path(dir_fig,"traj_seir2_paramB.pdf"),plot, width=40, height=40)

############
###########
SEIR2_I1=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("I1","I21"),to_sum="I1+I21",X="X_128.csv")
SEIR2_I2=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("I2","I12"),to_sum="I2+I12",X="X_128.csv")
SEIR2_I=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("I2","I12","I1","I21"),to_sum="I2+I12+I1+I21",X="X_128.csv")
SEIR2_S1=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("S1"),to_sum="S1/N",X="X_128.csv")
SEIR2_S2=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("S2"),to_sum="S2/N",X="X_128.csv")
SEIR2_S=proj(path1=pathSEIR2_psi, path2=pathSEIR2_psi_paramB,end_date = "2013-12-31",state=c("S"),to_sum="S/N",X="X_128.csv")

#Graphical settings
SEIR2_I1b=SEIR2_I1+ ylim(0,1000)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("A") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_I2b=SEIR2_I2+ ylim(0,1000)+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("B") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_Ib=SEIR2_I+ ylim(0,1000)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("C") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S1b=SEIR2_S1+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("D") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_S2b=SEIR2_S2+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("E") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))
SEIR2_Sb=SEIR2_S+ ylim(0,0.6)+theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle("F") +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))

#Combine plots
hlay <- rbind(c(1,2,3),
              c(4,5,6))

plot=grid.arrange(grobs=list(SEIR2_I1b,SEIR2_I2b,SEIR2_Ib,SEIR2_S1b,SEIR2_S2b,SEIR2_Sb), layout_matrix=hlay)
ggsave(file.path(dir_fig,"traj_seir2_psi_paramB.pdf"),plot, width=40, height=30)

