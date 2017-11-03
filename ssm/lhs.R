library(truncnorm)
library(rjson)
library(plyr)
library(dplyr)
library("lhs")

sample_from_lhs <- function(priors, theta_names=NULL,index) {
  
  theta <- rep(1, length(priors))
  names(theta) <- names(priors)
  prior_names <- names(priors)

  for(i in seq_along(theta)){
    
    prior <- priors[[prior_names[i]]]
    qprior <- paste0("q",prior$dist)
    prior_args <- c(prior$args,index[i])
    theta[i] <- do.call(qprior,prior_args)
  }
  
  return(theta)
}

generate_lhs= function(sample_size, priors){
  
  theta_estimated_names <- names(priors)

  LHS_square <- randomLHS(sample_size, length(theta_estimated_names))
  
  for(s in 1:sample_size){
    init_theta <- sample_from_lhs(priors, theta_estimated_names,index=LHS_square[s,])
    names(init_theta) <- theta_estimated_names
    init_theta_json =   list(name="values",description="initial values for the parameters",data=as.list(init_theta))
    
    covmat=diag(abs(init_theta)/10)
    colnames(covmat) = names(init_theta)
    rownames(covmat) = names(init_theta)
    covmat=as.list(as.data.frame(covmat))
    for(i in names(init_theta)){
        x <- covmat[[i]]
        x <- as.list(x)
        names(x) <- names(init_theta)
        covmat[[i]] <- x[x!=0] 
    }
    covmat_json=   list(name="covariance",description="covariance matrix",data=covmat)
    theta_json <- list(resources=list(init_theta_json, covmat_json))
    
    write(toJSON(theta_json),file=file.path(dir_lhs,paste0("theta_",s-1,".json")))
  }
}


do_lhs <- function(dir_lhs,dir_bin) {
  
  theta_files <- grep("theta_[0-9]+.json",list.files(dir_lhs),value=TRUE)
  
  for(i in seq_along(theta_files)){
    
    theta_file <- theta_files[i]
    cmd <- sprintf("cd %s; cat %s | ./simplex -M 10000 --prior > %s/theta_map_simplex_%s.json", 
                   dir_bin, 
                   file.path(paste0("../",dir_lhs),theta_file),
                   paste0("../",dir_lhs),
                   i-1
    )
    system(cmd, wait=FALSE)
    
  }
  
}

summarize_lhs <- function(dir_lhs, dir_model) {
  
  file_names <- grep(".*simplex.*",list.files(dir_lhs),value=TRUE)
  names(file_names) <- file_names
  
  df_ll = data.frame(matrix(NA,ncol = 2,nrow=length(file_names)))
  for(i in 1:length(file_names))  {
    res <- try(fromJSON(file=file.path(dir_lhs,file_names[i])), silent=FALSE)
    if(inherits(res, "try-error")){
      return(NA)
    } else {
      df_ll[i,]=c(file_names[i], res$resources[[3]]$data$log_ltp)		
    }
  }
  df_ll <- na.omit(df_ll)
  names(df_ll) <- c("file_name","log_ltp")
  df_ll$log_ltp=as.numeric(df_ll$log_ltp)
  
  df_theta <- ddply(df_ll,"file_name",function(df) {
    res <- fromJSON(file=file.path(dir_lhs,df$file_name))	
    unlist(res$resources[[1]]$data)
  }, .progress="text")
  
  df_lhs <- left_join(df_theta, df_ll, by="file_name")
  
  if(all(df_lhs$log_ltp <= 0)){
    df_lhs <- df_lhs %>% filter(log_ltp < -10)
  }

  df_lhs <- df_lhs %>% filter(log_ltp==max(log_ltp))
  
  cmd <- sprintf("cp %s %s/theta_map_simplex.json",file.path(dir_lhs, df_lhs$file_name),dir_model)
  system(cmd)
}



