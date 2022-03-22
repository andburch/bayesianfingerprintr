
getgamma <- function(meanz,sdz) {
  shape <- (meanz/sdz)^2
  rate<- meanz/sdz^2
  return(c(shape,rate))
  
}

extract_bayes <- function(samples, var) {
  returno <-NULL
  for (i in seq_along(samples)) {
    returno <- cbind(returno, samples[[i]][,var])
  }
  return(returno)
}

get_pseudopriors <- function (one.run, model.opt = c("2sex.2age","1sex.1age","2sex.1age","compare1v2age", "compare1v2ageNEXT")) {
  
  model1age<-NULL; model2age<-NULL;
  if (model.opt=="2sex.2age") {
    #turn mcmc.list into something easier
    mcmcMat= as.matrix(one.run)
    model2age$female$post_mean <- mcmcMat[,"hyper_age_mean[1]"]
    model2age$male$post_mean <- mcmcMat[,"hyper_age_mean[2]"]
    model2age$female$post_sd <- mcmcMat[,"hyper_age_sd[1]"]
    model2age$male$post_sd <- mcmcMat[,"hyper_age_sd[2]"]
    
    pseudopriorz <- 
      list(model2age = list(female = list(hyperage_mean2 = model2age$female$post_mean %>% {getgamma(mean(.),sd(.))},
                                          hyperage_sd2 = model2age$female$post_sd %>% {getgamma(mean(.),sd(.))}
      ),
      male =   list(hyperage_mean2 = model2age$male$post_mean %>% {getgamma(mean(.),sd(.))},
                    hyperage_sd2 = model2age$male$post_sd %>% {getgamma(mean(.),sd(.))}
      )
      )
      )
    
  } else if (model.opt == "2sex.1age") {
    #turn mcmc.list into something easier
    mcmcMat= as.matrix(one.run)
    
    model1age$post_mean <- mcmcMat[,"hyper_age_mean"]
    model1age$post_sd <- mcmcMat[,"hyper_age_sd"]
    pseudopriorz <- 
      list(model1age = list(hyperage_mean1 = model1age$post_mean %>% {getgamma(mean(.),sd(.))},
                            hyperage_sd1 = model1age$post_sd %>% {getgamma(mean(.),sd(.))}
      )
      )
    
    
  } else if (model.opt == "compare1v2age") {
    mcmcMat= as.matrix(one.run)
    
    m = mcmcMat[,"m"]
    model1age <- NULL; model2age <- NULL
    model1age$post_mean <- mcmcMat[,"hyper_age_mean1[1]"][m==1]
    model1age$post_sd <- mcmcMat[,"hyper_age_sd1[1]"][m==1]
    
    model2age$female$post_mean <- mcmcMat[,"hyper_age_mean2[1,2]"][m==2]
    model2age$male$post_mean <- mcmcMat[,"hyper_age_mean2[2,2]"][m==2]
    model2age$female$post_sd <- mcmcMat[,"hyper_age_sd2[1,2]"][m==2]
    model2age$male$post_sd <- mcmcMat[,"hyper_age_sd2[2,2]"][m==2]
    
    pseudopriorz <- 
      list(model1age = list(hyperage_mean1 = model1age$post_mean %>% {getgamma(mean(.),sd(.))},
                            hyperage_sd1 = model1age$post_sd %>% {getgamma(mean(.),sd(.))}
      ),
      model2age = list(female = list(hyperage_mean2 = model2age$female$post_mean %>% {getgamma(mean(.),sd(.))},
                                     hyperage_sd2 = model2age$female$post_sd %>% {getgamma(mean(.),sd(.))}
      ),
      male =   list(hyperage_mean2 = model2age$male$post_mean %>% {getgamma(mean(.),sd(.))},
                    hyperage_sd2 = model2age$male$post_sd %>% {getgamma(mean(.),sd(.))}
      )
      )
      )
  }
  return(pseudopriorz)
}




compare_1age2age <- function(db, iters=5000000, thin=100, n.chains=10, 
                             pre.iters=15000, pre.thin=3, pre.n.chains=3,
                             n.adapt=3000, n.burnin=8000, is.parallel=F) {
  
  if (n.burnin>=iters) stop("Burn-in must be less than the total iterations!")
  
  cat("Model with 2 sexes, 1 age distribution: \n")
  model1age <-  db %>%
    run_model("2sex.1age", "2sex.1age.txt", pre.iters, pre.thin,
              n.chains=pre.n.chains,return.raw.samples=T, is.parallel = is.parallel)
  
  cat("Model with 2 sexes, 2 age distributions: \n")
  model2age <- db %>%
    run_model("2sex.2age", "2sex.2age.txt", pre.iters, pre.thin,
              n.chains=pre.n.chains,return.raw.samples=T, is.parallel = is.parallel)
  
  pseudopriors <-
    c(get_pseudopriors(model1age, "2sex.1age"), get_pseudopriors(model2age, "2sex.2age"))
  
  model1age<-NULL; model2age<-NULL;
  
  
  saveRDS(pseudopriors, paste0(Sys.time() %>% str_sub(0,10),".pseudopriors.RDS"))
  
  second.run <- db %>%
    run_model("compare1v2ageNEXT", filename="compare1v2ageNEXT.txt",
              iters=iters,thin=thin,n.chains=10,return.raw.samples=T,pseudopriors=pseudopriors,
              n.burnin = n.burnin, n.adapt=n.adapt, is.parallel = is.parallel)
  
  effectiveSize(second.run)
  second.run %>% as.matrix() %>% .[,"m"] %>% table() %>% print() %>% cat()
  second.run %>% as.matrix() %>% .[,"m"] -> m_array
  
  #m_array is the vector of flagged values (1 or 2)
  prob_null_jags <- mean(m_array == 1)
  bf_jags <- prob_null_jags / (1 - prob_null_jags)
  cat("\n\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")
  cat(paste0("\n\nProb of NULL model (1 age distribution): ",prob_null_jags,".  Bayes Factor: ",bf_jags))
  
  return(second.run)
}



big.data <- read.csv("C:/Users/aburchil/Dropbox (ASU)/akiva-stuff/finaldata.csv", stringsAsFactors = F, na.strings=c("","NA")) %>% rowid_to_column()

noi <- 
  big.data %>% 
  mutate(groups = case_when(
    description=="ceramic" & phase %in% c("J") ~ "ceramic.J",
    description=="ceramic" & phase %in% c("H") ~ "ceramic.H",
    description=="ceramic" & phase %in% c("K","L") ~ "ceramic.LK",
    description=="mini"  ~ "mini",
    description=="figurine"  ~ "figurine",
    TRUE ~ description)
  ) %>% 
  mutate(mrb=lineup(mrb,.95)*mrb) %>% filter(groups=="ceramic.J")

for (i in unique(noi$groups)) {
  
  temp <- subset(noi, groups == i)
  if (nrow(temp) > 5) {
    
    cat(paste0("\n\n############",i,"#############\n#########################\n"))
    compare_1age2age(temp,50000,10)
    
  }
  
  
}





set.seed(137)
closer <- tribble(~diffz, ~modelz, ~probz, ~m.size)
temprun<-NULL

for (i in c(0,1,3,6)) {
  
  closer <-
    closer %>% 
    bind_rows(list(diffz = i,
                   modelz=list(simulate_data(100,.5,c(15+i,15-i),c(2,2), T)),
                   probz = NA,
                   m.size = NA
    )
    )
}

for (i in 1:nrow(closer)) {
  
  temprun <- compare_1age2age(closer$modelz[[i]], pre.iters=45000)
  m_array<- temprun %>% as.matrix() %>% .[,"m"] 
  
  #m_array is the vector of flagged values (1 or 2)
  closer$probz[[i]] <- mean(m_array == 1)
  closer$m.size[[i]] <- map(temprun, ~.[,"m"]) %>%
    as.mcmc.list() %>% effectiveSize() %>% as.numeric()
  
  temprun<-NULL
  
}








