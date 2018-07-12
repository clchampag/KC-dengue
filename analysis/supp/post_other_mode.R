library(gridExtra)  
require(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)

dir_fig = "./supp"
setwd("/ssm")


#setwd("/ssm")
ROOT = getwd() 
pathSEIR2 = paste(ROOT,"/SEIR2/mcmc", sep="")
pathSEIR2_paramB = paste(ROOT,"/SEIR2/mcmc-paramB", sep="")

pathSEIR2_psi = paste(ROOT,"/SEIR2-psi/mcmc", sep="")
pathSEIR2_psi_paramB = paste(ROOT,"/SEIR2-psi/mcmc-paramB", sep="")

N=161391

my_plot_posterior=function(model1,model2, param){
  
  plot_trace=read.csv(file = file.path(model1,"trace_1.csv"))[param]
  plot_trace2=read.csv(file = file.path(model2,"trace_1.csv"))[param]
  
  df_trace <- plot_trace %>% 
    gather_("theta","value",names(plot_trace)) 
  df_trace2 <- plot_trace2 %>% 
    gather_("theta","value",names(plot_trace2))

  p <- ggplot(df_trace, aes(fill=value))
  p <- p + geom_histogram(data=df_trace, aes(x=value, y=..density.. , fill="1"), 
                          alpha = 0.6 )
  p <- p + theme(strip.text.x = element_text(size = 55))
  p <- p + geom_histogram(data=df_trace2, aes(x=value, y=..density.. , fill="2"), 
                          alpha = 0.6 ) 
  p <- p + scale_fill_manual("",values=c("1"="dodgerblue","2"="red"), breaks=c("1","2"))
  p <- p + theme(panel.background = element_rect(fill = 'white', colour = 'white'))+ theme(axis.line = element_line(color="black", size = 1))
  p <- p + xlab("") + ylab("") + theme(legend.position="top") +theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))
  p <- p +  theme(legend.text=element_text(size=55), legend.key=element_rect(size=10,color="white"),legend.key.size=unit(2,"cm"))
  
  return(p)
} 


plot_fitness=my_plot_posterior(model1=pathSEIR2,model2=pathSEIR2_paramB,param="fitness")+scale_x_continuous(breaks=c(-2270,-2260,-2250))+
  ggtitle('1A') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+theme(legend.position="none")
  
plot_R0=my_plot_posterior(model1=pathSEIR2,model2=pathSEIR2_paramB,param="R0")+xlim(1, 3)+
  ggtitle('1B') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")

plot_rep=my_plot_posterior(model1=pathSEIR2,model2=pathSEIR2_paramB,param="rep_ndss")+scale_x_continuous(breaks=c(0.06,0.10,0.14), limits=c(0.05, 0.15))+
  ggtitle('1C') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")

#####
#SEIR2 psi
plot_psi_fitness=my_plot_posterior(model1=pathSEIR2_psi,model2=pathSEIR2_psi_paramB,param="fitness")+
  ggtitle('2A') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")

plot_psi_R0=my_plot_posterior(model1=pathSEIR2_psi,model2=pathSEIR2_psi_paramB,param="R0")+xlim(1, 3)+
  ggtitle('2B') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")

plot_psi_rep=my_plot_posterior(model1=pathSEIR2_psi,model2=pathSEIR2_psi_paramB,param="rep_ndss")+scale_x_continuous(breaks=c(0.06,0.10,0.14), limits=c(0.05, 0.15))+
  ggtitle('2C') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")

plot_psi_psi=my_plot_posterior(model1=pathSEIR2_psi,model2=pathSEIR2_psi_paramB,param="psi")+xlim(0.5, 0.9)+
  ggtitle('2D') +theme(plot.title=element_text( face="bold", size=65,hjust =0))+ theme(legend.position="none")



hlay <- rbind(c(1,2,3,NA),
              c(4,5,6,7))
seir2_post=grid.arrange(grobs=list(plot_fitness,plot_R0,plot_rep,
                                   plot_psi_fitness,plot_psi_R0,plot_psi_rep,plot_psi_psi),  layout_matrix=hlay)
ggsave(file.path(dir_fig,"post_paramB.pdf"),seir2_post, width=40, height=20, limitsize=F)
