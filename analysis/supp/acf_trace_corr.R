library(coda)
library(plyr)
library(reshape2)
library(lattice)
library(Cairo)
library(colorspace)
library(ggplot2)


dir_fig = "./supp"
setwd("/ssm")
ROOT = getwd() 

pathSEIR = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")

names_trace=function(mod, trace){
  if(mod=="pandey"){
    return(c("Hi(0)","Hs(0)","R0","\u03B2v","b","i","p","\u03A6","r","fitness", "index" ))
  }
  if(mod=="laneri"){
    return(c("Hi(0)","Hs(0)","R0","b","i","p","\u03A6","r","fitness", "index" ))
  }
  if(mod=="seir"){
    return(c("Hi(0)","Hs(0)","R0","b","i","p","\u03A6","r","fitness", "index" ))
  }
  if(mod=="seiar"){
    return(c("Hi(0)","Hs(0)","R0","b","i","p","\u03c1 a","\u03c1 h","\u03A6","fitness", "index" ))
  }
  if(mod=="seir2"){
    return(c("Hi1(0)","Hi2(0)","Hs(0)","Hs1(0)","Hs2(0)","R0","b","i","p","\u03A6","r","fitness", "index" ))
  }
  if(mod=="seir2psi"){
    return(c("Hi1(0)","Hi2(0)","Hs(0)","Hs1(0)","Hs2(0)","R0","b","i","p","\u03c8","\u03A6","r","fitness", "index" ))
  }
}


my_plot_trace=function(model, mcmc, name, mod){
  
  
  trace_files <- grep("trace_",list.files(paste(model,mcmc, sep="")),value=TRUE)
  trace=read.csv(file.path((paste(model,mcmc, sep="")),trace_files))
  names(trace)=names_trace(mod)
  trace_melt=melt(trace, id.vars = "index")
  trace_melt=trace_melt[trace_melt$variable !="fitness",]
  p= qplot(index, value, data=trace_melt, geom="line")+facet_grid(variable ~ ., scales="free")+
    theme(panel.background = element_rect(fill = 'white', colour = 'white'))+xlab("")+ylab("")+ 
    theme(strip.text.x = element_text(size = 40))
  p
}  

my_plot_correlation=function(model, mcmc, name, mod){
  trace_files <- grep("trace_",list.files(paste(model,mcmc, sep="")),value=TRUE)
  list_trace <- llply(file.path((paste(model,mcmc, sep="")),trace_files),read.csv, .progress="text")
  
  names(list_trace[[1]])=names_trace(mod)
  fitted_theta <- setdiff(names(list_trace[[1]]),c("fitness","index"))
  
    mycor=cor(list_trace[[1]][fitted_theta])
    mycor[upper.tri(mycor)]=NA
    diag(mycor)=0
    mycor=melt(mycor, na.rm = TRUE)
    ggplot(data = mycor, aes(Var2, Var1, fill = value))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, limit = c(-1,1), space = "Lab", name="") +
      theme_minimal()+xlab("")+ylab("")+ theme(axis.text = element_text(size=40), axis.text.x=element_text(angle=90))+
      theme(legend.text=element_text(size=25), legend.key=element_rect(size=10,color="white"),legend.key.size=unit(2,"cm"))
  
}  

#TRACE

CairoPNG(file.path(dir_fig,"trace_seir.png"),width = 480, height = 720)
my_plot_trace(model=pathSEIR,mcmc="/mcmc", mod="seir")+
  ggtitle('A') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPNG(file.path(dir_fig,"trace_laneri.png"),width = 480, height = 720)
my_plot_trace(model=pathLaneri,mcmc="/mcmc", mod="laneri")+
  ggtitle('B') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPNG(file.path(dir_fig,"trace_pandey.png"),width = 480, height = 720)
my_plot_trace(model=pathPandey,mcmc="/mcmc", mod="pandey")+
  ggtitle('C') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPNG(file.path(dir_fig,"trace_seiar.png"),width = 480, height = 720)
my_plot_trace(model=pathSEIAR,mcmc="/mcmc", mod="seiar")+
  ggtitle('D') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()


CairoPNG(file.path(dir_fig,"trace_seir2.png"),width = 480, height = 720)
my_plot_trace(model=pathSEIR2,mcmc="/mcmc", mod="seir2")+
  ggtitle('E') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()


CairoPNG(file.path(dir_fig,"trace_seir2psi.png"),width = 480, height = 720)
my_plot_trace(model=pathSEIR2_psi,mcmc="/mcmc", mod="seir2psi")+
  ggtitle('F') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

#CORRELATION
CairoPDF(file.path(dir_fig,"correlation_seir.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathSEIR,mcmc="/mcmc", mod="seir")+
  ggtitle('A') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPDF(file.path(dir_fig,"correlation_laneri.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathLaneri,mcmc="/mcmc", mod="laneri")+
  ggtitle('B') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPDF(file.path(dir_fig,"correlation_pandey.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathPandey,mcmc="/mcmc", mod="pandey")+
  ggtitle('C') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPDF(file.path(dir_fig,"correlation_seiar.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathSEIAR,mcmc="/mcmc", mod="seiar")+
  ggtitle('D') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPDF(file.path(dir_fig,"correlation_seir2.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathSEIR2,mcmc="/mcmc", mod="seir2")+
  ggtitle('E') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()

CairoPDF(file.path(dir_fig,"correlation_seir2psi.pdf"), width = 10, height = 10)
my_plot_correlation(model=pathSEIR2_psi,mcmc="/mcmc", mod="seir2psi")+
  ggtitle('F') +theme(plot.title=element_text( face="bold", size=65,hjust =0))
dev.off()
