library(coda)
library(plyr)
library(lattice)


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
    return(c("Hi(0)","Hs(0)","Vi(0)","R0","\u03B2v","b","i","p","r","fitness", "index" ))
  }
  if(mod=="laneri"){
    return(c("Hi(0)","L(0)","Hs(0)","R0","b","i","p","r","fitness", "index" ))
  }
  if(mod=="seir"){
    return(c("Hi(0)","Hs(0)","R0","b","i","p","r","fitness", "index" ))
  }
  if(mod=="seiar"){
    return(c("Hi(0)","Hs(0)","R0","b","i","p","\u03c1 a","\u03c1 h","fitness", "index" ))
  }
  if(mod=="seir2"){
    return(c("Hi1(0)","Hi2(0)","Hs(0)","Hs1(0)","Hs2(0)","R0","b","i","p","r","fitness", "index" ))
  }
  if(mod=="seir2psi"){
    return(c("Hi1(0)","Hi2(0)","Hs(0)","Hs1(0)","Hs2(0)","R0","b","i","p","\u03c8","r","fitness", "index" ))
  }
}


my_plot_trace=function(model, mcmc, name, mod){
  trace_files <- grep("trace_",list.files(paste(model,mcmc, sep="")),value=TRUE)
  list_trace <- llply(file.path((paste(model,mcmc, sep="")),trace_files),read.csv, .progress="text")
  length_trace <- sapply(list_trace, nrow)
  add_to_index <- c(0,cumsum(length_trace)[-length(length_trace)])
  for(chain in seq_along(add_to_index)){
    list_trace[[chain]]$index <- list_trace[[chain]]$index + add_to_index[chain]
  }	
  
  names(list_trace[[1]])=names_trace(mod)
  fitted_theta <- setdiff(names(list_trace[[1]]),c("fitness","index"))
  my_mcmc <- mcmc.list(llply(list_trace,function(x) {mcmc(x[fitted_theta])}))
  
  p <- xyplot(my_mcmc,strip = FALSE, strip.left = strip.custom(style = 1, horizontal = FALSE))
  p
  
}  


#TRACE
my_plot_trace(model=pathSEIR,mcmc="/mcmc", mod="seir")
my_plot_trace(model=pathLaneri,mcmc="/mcmc", mod="laneri")
my_plot_trace(model=pathPandey,mcmc="/mcmc", mod="pandey")
my_plot_trace(model=pathSEIAR,mcmc="/mcmc", mod="seiar")
my_plot_trace(model=pathSEIR2,mcmc="/mcmc", mod="seir2")
my_plot_trace(model=pathSEIR2_psi,mcmc="/mcmc", mod="seir2psi")
