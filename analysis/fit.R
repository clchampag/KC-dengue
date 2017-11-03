require(data.table)
library(ggplot2)
library(gridExtra)  



setwd("/ssm")
ROOT = getwd() 
pathSEIR1 = paste(ROOT,"/SEIR", sep="")
pathLaneri = paste(ROOT,"/Laneri", sep="")
pathPandey = paste(ROOT,"/Pandey", sep="")
pathSEIAR = paste(ROOT,"/SEIAR", sep="")
pathSEIR2 = paste(ROOT,"/SEIR2", sep="")
pathSEIR2_psi = paste(ROOT,"/SEIR2-psi", sep="")

fit_model=function(path, incidence, start_date,end_date, file_data=NULL){
  fileNames = list.files(file.path(path,"mcmc/"), pattern = "X_128", full.names = T)
  fileDf = fread(fileNames[1], sep = ",", select = c("date", paste0("ran_", incidence),"index"))
  names(fileDf)=c("date", "ran_incidence","index")
  
  fileData=list.files(file.path(path,"data"), pattern = "data.csv", full.names = T)
  if (is.null(file_data)==FALSE){fileData=list.files(datapath, pattern = "data_proj.csv", full.names = T)}
  data <- fread(fileData, sep = ",", select = c("date", incidence))
  names(data)=c("date", "data")
  data=as.data.frame(data)
  data$date=as.Date(data$date)
  
  median=fileDf[, median(ran_incidence), by = c("date")]
  q975=fileDf[, quantile(ran_incidence, probs=c(0.975)), by = c("date")]
  q025=fileDf[, quantile(ran_incidence, probs=c(0.025)), by = c("date")]
  
  qmin=fileDf[, quantile(ran_incidence, probs=c(0.001)), by = c("date")]
  qmax=fileDf[, quantile(ran_incidence, probs=c(0.999)), by = c("date")]
  
  tot=data.frame(date=as.Date(median$date), median=median$V1, q975=q975$V1, q025=q025$V1)
  tot=tot[as.Date(tot$date)<end_date & as.Date(tot$date)> start_date,]
  tot=merge(data,tot, all=T)
  tot=tot[as.Date(tot$date)<end_date & as.Date(tot$date)> start_date,]
  
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  ggplot(tot,aes_string(x='date',y ="median")) +
    geom_ribbon(aes_string(ymin= "q025", ymax="q975"),fill="skyblue",alpha=0.55) +
    #geom_ribbon(aes_string(ymin= "lower_50", ymax="upper_50"),fill="skyblue",alpha=1) +
    geom_line(color="dodgerblue4",size=1)+ theme(axis.line.x = element_line(color="black"), axis.line.y = element_line(color="black"))+
    geom_point(aes_string(x='date', y="data"),color="black", size=5)+
    theme(axis.text.x = element_text(size=45)) + theme(axis.text.y = element_text(size=45))+theme(strip.background = element_blank())+
    theme(axis.title = element_text(size=16)) +ylab("") + xlab("") +
    theme(panel.background = element_rect(fill = 'white'), panel.grid.major = element_line(colour = "lightgrey", linetype = "dashed"))
}
################
#FIT ON NDSS DATA
# Plot each model separately
SEIR=fit_model(pathSEIR,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")
Laneri=fit_model(pathLaneri,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")
Pandey=fit_model(pathPandey,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")
SEIR2=fit_model(pathSEIR2,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")
SEIR2psi=fit_model(pathSEIR2_psi,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")
SEIAR=fit_model(pathSEIAR,incidence="incidence_ndss_obs",start_date="2002-01-01",end_date="2013-12-31")

# Graphical settings
SEIR=SEIR+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 
Laneri=Laneri+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Laneri') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 
Pandey=Pandey+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Pandey') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 
SEIAR=SEIAR+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIAR') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 
SEIR2=SEIR2+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 
SEIR2psi=SEIR2psi+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2 psi') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,250) 

#Compbine plots
hlay <- rbind(c(1,2,3),
              c(4,5,6))

plot_all=grid.arrange(grobs=list(SEIR,Laneri,Pandey,SEIAR,
                                 SEIR2, SEIR2psi),  layout_matrix=hlay)


################
#FIT ON DENFREE DATA
# Plot each model separately
SEIR=fit_model(pathSEIR1,incidence="incidence_denfree_obs",start_date = "2012-01-01",end_date="2013-12-31")
Laneri=fit_model(pathLaneri,incidence="incidence_denfree_obs",start_date = "2012-01-01",end_date="2013-12-31")
Pandey=fit_model(pathPandey,incidence="incidence_denfree_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIAR_s=fit_model(pathSEIAR,incidence="incidence_denfree_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIAR_a=fit_model(pathSEIAR,incidence="incidence_asympto_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIR2_1=fit_model(pathSEIR2,incidence="incidence_D1_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIR2_2=fit_model(pathSEIR2,incidence="incidence_D2_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIR2psi_1=fit_model(pathSEIR2_psi,incidence="incidence_D1_obs",start_date = "2012-01-01",end_date="2013-12-31")
SEIR2psi_2=fit_model(pathSEIR2_psi,incidence="incidence_D2_obs",start_date = "2012-01-01",end_date="2013-12-31")

# Graphical settings
SEIR=SEIR+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90)  +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
Laneri=Laneri+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Laneri') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90)  +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
Pandey=Pandey+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Pandey') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90)  +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIR2_1=SEIR2_1+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2: D1') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90) +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIR2_2=SEIR2_2+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2: D2') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ylim(0,40)+scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIAR_s=SEIAR_s+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIAR: sympt.') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90) +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIAR_a=SEIAR_a+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIAR: asympt.') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIR2psi_1=SEIR2psi_1+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2 psi: D1') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,90) +scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")
SEIR2psi_2=SEIR2psi_2+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2 psi: D2') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ylim(0,40)+scale_x_date(date_breaks=("1 year") ,date_labels = "%Y")

#Combine plots
hlay <- rbind(c(1,2,3),
              c(4,5,NA),
              c(6,7,NA),
              c(8,9,NA))

plot_all=grid.arrange(grobs=list(SEIR,Laneri,Pandey,SEIAR_s,SEIAR_a,
                                 SEIR2_1,SEIR2_2,SEIR2psi_1,SEIR2psi_2),  layout_matrix=hlay)



#############
# Fit on projections (2014-2015)
datapath = "../analysis - copie/"

# Plot each model separately
SEIR=fit_model(pathSEIR,incidence="incidence_ndss_obs", end_date = "2015-12-31",start_date = "2014-01-01" ,file_data = datapath)
Laneri=fit_model(pathLaneri,incidence="incidence_ndss_obs",start_date="2014-01-01",end_date="2015-12-31", file_data = datapath)
Pandey=fit_model(pathPandey,incidence="incidence_ndss_obs", end_date = "2015-12-31",start_date = "2014-01-01" , file_data = datapath)
SEIR2=fit_model(pathSEIR2,incidence="incidence_ndss_obs",end_date = "2015-12-31",start_date = "2014-01-01", file_data = datapath)
SEIAR=fit_model(pathSEIAR,incidence="incidence_ndss_obs",end_date = "2015-12-31",start_date = "2014-01-01" , file_data = datapath)
SEIR2_psi=fit_model(pathSEIR2_psi,incidence="incidence_ndss_obs", end_date = "2015-12-31",start_date = "2014-01-01", file_data = datapath)



#Graphical settings
SEIR=SEIR+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 
Laneri=Laneri+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Laneri') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 
Pandey=Pandey+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('Pandey') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 
SEIR2=SEIR2+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 
SEIR2_psi=SEIR2_psi+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIR2 psi') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 
SEIAR=SEIAR+ theme(plot.margin=unit(c(1,1,0.5,0.5), "cm"))+ ggtitle('SEIAR') +theme(plot.title=element_text( face="bold", size=80,hjust =0.1))+ ylim(0,110) 

#Combine plots
hlay <- rbind(c(1,2,3),
              c(4,5,6))

plot_all=grid.arrange(grobs=list(SEIR,Laneri,Pandey,SEIAR,
                                 SEIR2, SEIR2_psi),  layout_matrix=hlay)

