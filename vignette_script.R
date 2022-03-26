#' ---
#' title: "Bayesian Mixture-Model Analysis of Ancient Fingerprints"
#' author: "Andrew Burchill"
#' date: "2020-22-02"
#' output:
#'    html_document:
#'       self_contained: true
#'       theme: lumen
#'       includes:
#'         after_body: footer.html
#'         before_body: style.html
#'       code_folding: "hide"
#'       toc: true
#'       toc_float: 
#'         collapsed: true
#'         smooth_scroll: true
#'       
#' always_allow_html: yes  
#' ---
#' 

#+ setup, include=FALSE, results='hide'

require("dplyr")
require("tibble")
require("rjags")
require("magrittr")
require("stringr")
require("purrr")
require("ggplot2")
require("tidyr")
require("jagsUI")


fingerprint_file_path <- "kralik_novotny_fingerprints.csv"


lineup <- function(mrbs,x) {
  ref_data <-
    read.csv(file=fingerprint_file_path) %>%
    as_tibble() %>% 
    group_by(ID) %>%
    summarize(mrb=mean(RB), age=unique(age), sex=unique(sex))
  
  quantile(ref_data$mrb,x)/quantile(mrbs,x)
}

lineup2 <- function(mrbs,x,type=c("males","adult males", "adults"),adult_age=15) {
  ref_data<-    read.csv(file=fingerprint_file_path) %>%
    as_tibble()
  
  ref_data <-
    ref_data %>% 
    group_by(ID) %>%
    summarize(mrb=mean(RB), age=unique(age), sex=unique(sex)) 
  db = switch(type,
              "adults" = filter(ref_data, age >= adult_age),
              "males" = filter(ref_data, sex == 2),
              "adult males" = filter(ref_data, age >= adult_age, sex == 2)
  )
  
  quantile(db$mrb,x)/quantile(mrbs,x)
}


make_modelstring <-
  function(model.opt = c("2sex.2age","1sex.1age","2sex.1age","1sex.2age", "compare1v2age", "compare1v2ageNEXT"),
           filename="bayesmodel.txt", pseudopriors=NULL) {
    
    raw.string.head<-NULL
    raw.string.tail <-NULL
    raw.string <-NULL
    single.sex.opts <- c("1sex.1age")
    ## Men are SEX = 2!!! Women are SEX = 1
    {  
      raw.string.head <-
        "model {
    
#stdev for likelihood function
sd ~ dgamma(2, 10)

for (i in 1:N) {

  #Likelihood function is a normal distribution
  mrb[i] ~ dnorm(u[i], 1/sd^2)
  
  #calculate log... likelihood?
  LogLik[i] <- logdensity.norm(mrb[i],u[i], 1/sd^2)
  
  #If you want the MRB likelihood function to be gamma distributed, use this:
  #mrb[i] ~ dgamma(shape[i], rate[i])
  #shape[i] <- (u[i]/sd)^2
  #rate[i] <- u[i]/sd^2
  
  #The asymptote equation is: Asym+(R0-Asym)*exp(-exp(lrc)*input)
   "
      
      raw.string.tail <- NULL
      nodes <- NULL
      
      {
        nodes["2sex.2age"] <-
          list(c("age","sex","mrb","sd","LogLik",
                 "lrc[1]","Asym[1]","R0[1]","lrc[2]","Asym[2]","R0[2]","pClust",
                 "hyper_age_mean[2]","hyper_age_mean[1]","hyper_age_sd[2]","hyper_age_sd[1]"
          ))
        
        raw.string.tail["2sex.2age"] <-
          "u[i] <- (Asym[sex[i]]) +
          (R0[sex[i]] - (Asym[sex[i]]))*
          exp(-1*exp(lrc[sex[i]]) * age[i])
  
  #The latent 'age' variable is pulled from one of two gamma distributions,
  #depending on the inferred sex. The gamma distributions for each sex have their
  #own hyper priors, allowing them to be updated based on the data.
  
  age[i] ~ dgamma((hyper_age_mean[sex[i]]/hyper_age_sd[sex[i]])^2, hyper_age_mean[sex[i]]/hyper_age_sd[sex[i]]^2)
  
  #Sex is chosen based off the probability of belonging to one or another distribution (sex)
  sex[i] ~ dcat(sexprobs)
  
}

#hyperpriors for the gamma distribution of ages for each sex
hyper_age_mean[1] ~ dgamma(4,.2)
hyper_age_mean[2] ~ dgamma(4,.2)

hyper_age_sd[1] ~ dgamma(1,.1)
hyper_age_sd[2] ~ dgamma(1,.1)


#Normal, but TIGHT priors for the asymptote parameters,
#derived from using an asymptotic regression model on Kralik's data.
Asym[2] ~ dnorm(%s,1/.008^2)
R0[2] ~ dnorm(%s,1/.008^2)
lrc[2] ~ dnorm(%s,1/.03^2)

Asym[1] ~ dnorm(%s,1/.008^2)
R0[1] ~ dnorm(%s,1/.008^2)
lrc[1] ~ dnorm(%s,1/.03^2)


#the sex distribution
pClust ~ dbeta(1,1)
sexprobs[1] = pClust #female
sexprobs[2] = 1-pClust #male

}
"
      }

{
  nodes["1sex.2age"] <-
    list(c("age","sex","mrb","sd","LogLik",
           "lrc","Asym","R0","pClust",
           "hyper_age_mean[2]","hyper_age_mean[1]","hyper_age_sd[2]","hyper_age_sd[1]"
    ))
  
  raw.string.tail["1sex.2age"] <-
    "u[i] <- (Asym) +
          (R0 - (Asym))*
          exp(-1*exp(lrc) * age[i])
  
  #The latent 'age' variable is pulled from one of two gamma distributions,
  #depending on the inferred age group (referred to as `sex` for consistency 
  #with other models) The gamma distributions for each `sex` have their
  #own hyper priors, allowing them to be updated based on the data.
  
  age[i] ~ dgamma((hyper_age_mean[sex[i]]/hyper_age_sd[sex[i]])^2, hyper_age_mean[sex[i]]/hyper_age_sd[sex[i]]^2)
  
  #`Sex` is chosen based off the probability of belonging to one or another distribution (young v old)
  sex[i] ~ dcat(sexprobs)
  
}

#hyperpriors for the gamma distribution of ages for each `sex`
hyper_age_mean[1] ~ dgamma(4,.2)

#also, the mean age difference should maybe be AROUND ~3-13 years?
diffo ~ dgamma(3.2,0.4)
hyper_age_mean[2] <- hyper_age_mean[1] + diffo

hyper_age_sd[1] ~ dgamma(1,.1)
hyper_age_sd[2] ~ dgamma(1,.1)


#Normal, but TIGHT priors for the asymptote parameters,
#derived from using an asymptotic regression model on Kralik's data.
Asym ~ dnorm(%s,1/.008^2)
R0 ~ dnorm(%s,1/.008^2)
lrc ~ dnorm(%s,1/.03^2)


#the sex distribution
pClust ~ dbeta(1,1)
sexprobs[1] = pClust #female
sexprobs[2] = 1-pClust #male

}
"
}

{
  nodes["1sex.1age"] <-
    list(c("age","mrb","sd","LogLik",
           "lrc","Asym","R0",
           "hyper_age_mean","hyper_age_sd"
    ) )
  
  raw.string.tail["1sex.1age"] <-
    "u[i] <- (Asym) +
          (R0 - (Asym))*
          exp(-1*exp(lrc) * age[i])
  
  #The latent 'age' variable is pulled from a gamma distributions with
  #hyper priors, allowing it to be updated based on the data.
  
  age[i] ~ dgamma((hyper_age_mean/hyper_age_sd)^2, hyper_age_mean/hyper_age_sd^2)
  
  
}

#hyperpriors for the gamma distribution of ages for males
hyper_age_mean ~ dgamma(4,.2)

hyper_age_sd ~ dgamma(1,.1)


#Normal, but TIGHT priors for the asymptote parameters,
#derived from using an asymptotic regression model on Kralik's data.

Asym ~ dnorm(%s,1/.008^2)
R0 ~ dnorm(%s,1/.008^2)
lrc ~ dnorm(%s,1/.03^2)


}
"
}

{
  nodes["2sex.1age"] <-
    list(c("age","sex","mrb","sd","LogLik",
           "lrc[1]","Asym[1]","R0[1]","lrc[2]","Asym[2]","R0[2]","pClust",
           "hyper_age_mean","hyper_age_sd"
    ))
  
  raw.string.tail["2sex.1age"] <-
    "u[i] <- (Asym[sex[i]]) +
          (R0[sex[i]] - (Asym[sex[i]]))*
          exp(-1*exp(lrc[sex[i]]) * age[i])
  
  #The latent 'age' variable is pulled from a gamma distribution with its
  #own hyper priors, allowing it to be updated based on the data.
  
 age[i] ~ dgamma((hyper_age_mean/hyper_age_sd)^2, hyper_age_mean/hyper_age_sd^2)
  
  
  #Sex is chosen based off the probability of belonging to one or another distribution (sex)
  sex[i] ~ dcat(sexprobs)
  
}

#hyperpriors for the gamma distribution of ages 
hyper_age_mean ~ dgamma(4,.2)

hyper_age_sd ~ dgamma(1,.1)

#Normal, but TIGHT priors for the asymptote parameters,
#derived from using an asymptotic regression model on Kralik's data.
Asym[2] ~ dnorm(%s,1/.008^2)
R0[2] ~ dnorm(%s,1/.008^2)
lrc[2] ~ dnorm(%s,1/.03^2)

Asym[1] ~ dnorm(%s,1/.008^2)
R0[1] ~ dnorm(%s,1/.008^2)
lrc[1] ~ dnorm(%s,1/.03^2)


#the sex distribution
pClust ~ dbeta(1,1)
sexprobs[1] = pClust #female
sexprobs[2] = 1-pClust #male

}
"
}

{
  nodes["compare1v2age"] <-
    list(c("age","sex","mrb","sd","LogLik",
           "lrc[1]","Asym[1]","R0[1]","lrc[2]","Asym[2]","R0[2]","pClust",
           "hyper_shape1","hyper_rate1", "m",
           "hyper_shape2[1]","hyper_rate2[1]","hyper_shape2[2]","hyper_rate2[2]"
    ))
  
  raw.string.tail["compare1v2age"] <-
    "u[i] <- (Asym[sex[i]]) +
          (R0[sex[i]] - (Asym[sex[i]]))*
          exp(-1*exp(lrc[sex[i]]) * age[i])
  
  #The latent 'age' variable is pulled from one of two gamma distributions,
  #depending on the inferred sex. The gamma distributions for each sex have their
  #own hyper priors, allowing them to be updated based on the data.
  
  age[i] ~ dgamma(age_shape[i], age_rate[i])
  
  age_shape[i] = equals(m,1)*hyper_shape1 +  equals(m,2)*hyper_shape2[sex[i]]
  age_rate[i] =  equals(m,1)*hyper_rate1 +  equals(m,2)*hyper_rate2[sex[i]]
  
  #Sex is chosen based off the probability of belonging to one or another distribution (sex)
  sex[i] ~ dcat(sexprobs)
  
}


#Normal, but TIGHT priors for the asymptote parameters,
#derived from using an asymptotic regression model on Kralik's data.


Asym[2] ~ dnorm(%s,1/.008^2)
R0[2] ~ dnorm(%s,1/.008^2)
lrc[2] ~ dnorm(%s,1/.03^2)

Asym[1] ~ dnorm(%s,1/.008^2)
R0[1] ~ dnorm(%s,1/.008^2)
lrc[1] ~ dnorm(%s,1/.03^2)

    
    
#which model is 'active'?
#m=1 -> 1 age dist, m=2 -> 2 age dists
m ~ dcat(mPriorProb[]) 

mPriorProb[1] <- .5
mPriorProb[2] <- .5


#hyper age 1 is the female+male age distribution
hyper_age_mean1_prior ~ dgamma(4,.2)
hyper_age_sd1_prior ~ dgamma(1,.1)
#pseudopriors
hyper_age_mean1_pseudoprior ~ dgamma(%s,%s)
hyper_age_sd1_pseudoprior ~ dgamma(%s,%s)

hyper_age_mean1 = equals(m,1)*hyper_age_mean1_prior + equals(m,2)*hyper_age_mean1_pseudoprior
hyper_age_sd1 = equals(m,1)*hyper_age_sd1_prior + equals(m,2)*hyper_age_sd1_pseudoprior
hyper_shape1 = (hyper_age_mean1/hyper_age_sd1)^2
hyper_rate1  = hyper_age_mean1/(hyper_age_sd1)^2


      
#hyperpriors for the gamma distribution of ages for 2 sexes
#2=male, 1=female
#old format [sex, model]
hyper_age_mean2_prior[1] ~ dgamma(4,.2)
hyper_age_mean2_prior[2] ~ dgamma(4,.2)
hyper_age_sd2_prior[1] ~ dgamma(1,.1)
hyper_age_sd2_prior[2] ~ dgamma(1,.1)
#pseudopriors
hyper_age_mean2_pseudoprior[1] ~ dgamma(%s,%s)
hyper_age_mean2_pseudoprior[2] ~ dgamma(%s,%s)
hyper_age_sd2_pseudoprior[1] ~ dgamma(%s,%s)
hyper_age_sd2_pseudoprior[2] ~ dgamma(%s,%s)



hyper_age_mean2[1] = equals(m,2)*hyper_age_mean2_prior[1] + equals(m,1)*hyper_age_mean2_pseudoprior[1]
hyper_age_sd2[1] = equals(m,2)*hyper_age_sd2_prior[1] + equals(m,1)*hyper_age_sd2_pseudoprior[1]
hyper_shape2[1] = (hyper_age_mean2[1]/hyper_age_sd2[1])^2
hyper_rate2[1]  = hyper_age_mean2[1]/(hyper_age_sd2[1])^2

hyper_age_mean2[2] = equals(m,2)*hyper_age_mean2_prior[2] + equals(m,1)*hyper_age_mean2_pseudoprior[2]
hyper_age_sd2[2] = equals(m,2)*hyper_age_sd2_prior[2] + equals(m,1)*hyper_age_sd2_pseudoprior[2]
hyper_shape2[2] = (hyper_age_mean2[2]/hyper_age_sd2[2])^2
hyper_rate2[2]  = hyper_age_mean2[2]/(hyper_age_sd2[2])^2


#the sex distribution
pClust ~ dbeta(1,1)
sexprobs[1] = pClust #female
sexprobs[2] = 1-pClust #male

}
"
}

#if you're updating pseudopriors THEN do this
if (!is.null(pseudopriors)){
  nodes["compare1v2ageNEXT"] <- nodes["compare1v2age"];
  raw.string.tail["compare1v2ageNEXT"]  <-raw.string.tail["compare1v2age"]
} 

    }

#this is a list of models, without the right coefficients inserted
raw.string <- paste(raw.string.head, raw.string.tail) %>% 
  setNames(names(raw.string.tail)) 

## Men are SEX = 2!!! Women are SEX = 1
ref_data<-read.csv(file=fingerprint_file_path) %>% as_tibble()

ref_data <-
  ref_data %>% 
  group_by(ID) %>%
  summarize(mrb=mean(RB), age=unique(age), sex=unique(sex)) 

male_coef <-
  nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(ref_data, sex==2)) %>% 
  coef()
female_coef <-   
  nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(ref_data, sex==1)) %>% 
  coef()


insert_coefficients <- function(stringz, namez) {
  
  if (namez %in% single.sex.opts) {
    sprintf(stringz, 
            male_coef['Asym'],
            male_coef['r0'],
            male_coef['lrc']
    ) %>% return()
    
  } else if (namez %in% c("compare1v2age")) {
    sprintf(stringz, 
            male_coef['Asym'],
            male_coef['r0'],
            male_coef['lrc'],
            female_coef['Asym'],
            female_coef['r0'],
            female_coef['lrc'],
            4,.2,1,.1,4,.2,4,.2,1,.1,1.,1  #all the NORMAL priors
    ) %>% return()
    
  } else if(namez %in% c("compare1v2ageNEXT")) {
    sprintf(stringz, 
            male_coef['Asym'],
            male_coef['r0'],
            male_coef['lrc'],
            female_coef['Asym'],
            female_coef['r0'],
            female_coef['lrc'],
            pseudopriors$model1age$hyperage_mean1[1],
            pseudopriors$model1age$hyperage_mean1[2],
            pseudopriors$model1age$hyperage_sd1[1],
            pseudopriors$model1age$hyperage_sd1[2],
            pseudopriors$model2age$female$hyperage_mean2[1],
            pseudopriors$model2age$female$hyperage_mean2[2],
            pseudopriors$model2age$male$hyperage_mean2[1],
            pseudopriors$model2age$male$hyperage_mean2[2],
            pseudopriors$model2age$female$hyperage_sd2[1],
            pseudopriors$model2age$female$hyperage_sd2[2],
            pseudopriors$model2age$male$hyperage_sd2[1],
            pseudopriors$model2age$male$hyperage_sd2[2]
    ) %>% return()
  } else {
    sprintf(stringz, 
            male_coef['Asym'],
            male_coef['r0'],
            male_coef['lrc'],
            female_coef['Asym'],
            female_coef['r0'],
            female_coef['lrc']
    ) %>% return()
    
  }                              
  
}
#cat(names(raw.string))

#makes final versions of all models, with the right coefficients
model.string <- purrr::imap_chr(raw.string, insert_coefficients)

#this saves the model as an external text file (useful for multi-cores)
sink(filename)
cat(model.string[model.opt],
    fill = TRUE)
sink()

#return the parameters you need
return(nodes[model.opt])
  }


run_model <- function(  
  db, model.opt = c("2sex.2age","1sex.1age","2sex.1age","1sex.2age", "compare1v2age", "compare1v2ageNEXT"), filename="bayesmodel.txt",
  iters, thin, is.parallel=F, DIC=F, n.chains=10, return.raw.samples=F,
  n.adapt=3000, n.burnin=100, fix.biggest=T,pseudopriors=NULL){
  
  single.sex.model <- model.opt %in% c("1sex.1age") #also found in make_modelstring
  
  nodes <- make_modelstring(model.opt, filename, pseudopriors)
  
#add a rowid. column to make analyzing the results easier
  if("rowid." %in% colnames(db)) {
    stop("Supplied dataframe should not have a 'rowid.' column!")
  } else {
    db <- db %>% tibble::rowid_to_column("rowid.")
  }
#give the model run a name, hopefully based off a "groups" column
  
  if("groups" %in% colnames(db)) {
    run.name <- unique(db$groups)
  } else {
    run.name <- "bayes_run"
  }
  
  #make data ready for rjags
  temp.data <-
    db %>% 
    transmute(mrb=mrb, sex=NA) %>%
    as.list() %>% {`[<-`(., "N", length(.$mrb))} 
  
  #Set the largest fingerprint size to be male (this is an assumption)
  if(fix.biggest){
    temp.data$sex[which.max(temp.data$mrb)]=2 #2= male, 1=female
  } else {
    temp.data$sex<-NULL
  }
  
  #UNLESS THERE'S ONLY ONE GROUP
  if (model.opt=="1sex.1age") {
    temp.data$sex <-NULL
  }
  
  #RUN THE JAGS PART
  print(now <-Sys.time())
  #print(temp.data$N)
  
  #determine if you want multiple cores or not
  
  if (is.parallel) {
    samps <- jagsUI::jags(data=temp.data, inits=NULL, parameters.to.save = nodes %>% purrr::flatten_chr(), 
                          model.file = filename, n.chains=n.chains, n.cores=n.chains, DIC=DIC,
                          parallel = is.parallel, n.adapt=n.adapt, n.burnin=n.burnin,
                          n.thin=thin, n.iter = iters)
  } else {
    samps <- jagsUI::jags(data=temp.data, inits=NULL, parameters.to.save = nodes %>% purrr::flatten_chr(), 
                          model.file = filename, n.chains=n.chains, n.cores=NULL, DIC=DIC,
                          parallel = is.parallel, n.adapt=n.adapt, n.burnin=n.burnin,
                          n.thin=thin, n.iter = iters)
  }
  
  #remove all the other stuff jagsUI returns (wish I had used jagsUI earlier...)
  samps <- samps$samples
  print(Sys.time()-now)
  attr(samps,"original.group") <- run.name
  
  
  #return raw samples? Or return tidy posterior dataframe?
  
  if (return.raw.samples == T) {
    return(samps)
  } else {
    #title is either based on the "orgiginal group" or is named "bayes"
    title <- ifelse(is.null(attributes(samps)$original.group),"bayes",attributes(samps)$original.group)
    
    post_data <- 
      tibble(est_age = purrr::map_df(samps, ~as_tibble(.)  %>% 
                                select(starts_with("age")), simplify=T) %>% unlist(),
             
             LogLik = purrr::map_df(samps, ~as_tibble(.)  %>% 
                               select(starts_with("LogLik")), simplify=T) %>% unlist(),
             
             #List estimate sex as only males NOW, but then we'll fix it later in the pipe
             #it had some weird behavior
             est_sex = purrr::map_df(samps, ~as_tibble(.)  %>% 
                                select(starts_with("age")), simplify=T) %>%
               unlist() %>% length() %>% rep(2,.),
             
             rowid. = purrr::map_df(samps, ~as_tibble(.)  %>% 
                               select(starts_with("age")), simplify=T) %>%
               unlist() %>% 
               names() %>% 
               stringr::str_replace_all("(age[:punct:])([:digit:]+)([:punct:][:digit:]+)","\\2") %>% 
               as.numeric(),
             type = title
      ) %>% 
      
      #add actual sex estimates when needed
      
      {if(!single.sex.model) mutate(.,
                                    est_sex = purrr::map_df(samps, ~as_tibble(.) %>% select(starts_with("sex")),
                                                     simplify=T) %>% unlist()
      ) else .} %>% 
      
      mutate(est_sex = case_when(est_sex == 2 ~ "male", est_sex == 1 ~ "female", 
                                 TRUE ~ as.character(est_sex)) %>% factor(levels=c("male","female"))
      ) %>% 
      left_join(db, by = "rowid.")
    attr(post_data,"original.group") <- run.name
    return(post_data)
  }
}


run_single_model <- function(db, iters=1500, thin=3, return.raw.samples=F) {
  
  run_model(db, model.opt = "2sex.2age", filename="2sex.2age.txt",
            iters=iters, thin=thin, n.chains=3, return.raw.samples = return.raw.samples)
  
}


run_1age_model <- function(db, iters=1500, thin=3, return.raw.samples=F) {
  
  run_model(db, model.opt = "2sex.1age", filename="2sex.1age.txt",
            iters=iters, thin=thin, n.chains=3, return.raw.samples = return.raw.samples)
  
  
  
  
}


run_MEN_ONLY_single_model <- function(db, iters=1500, thin=3, return.raw.samples=F) {
  
  run_model(db, model.opt = "1sex.1age", filename="1sex.1age.txt",
            iters=iters, thin=thin, n.chains=3, return.raw.samples = return.raw.samples)
  
  
}




posterior_plot <- function(db, is.it.raw.samples=F, percentilex=NULL) {
  
  require("ggplot2")
  title <- ifelse(is.null(attributes(db)$original.group),"bayes",attributes(db)$original.group)
  
  post_data <- db
  
  
  if (is.it.raw.samples==T){
    
    post_data <- NULL
    post_data <- 
      tibble(est_age = purrr::map_df(db, ~as_tibble(.)  %>% 
                                select(starts_with("age")), simplify=T) %>% unlist(),
             est_sex = purrr::map_df(db, ~as_tibble(.)  %>% 
                                select(starts_with("sex")), simplify=T) %>% unlist()
      ) 
    
  } 
  
  p <- 
    ggplot(post_data, aes(x=est_age, color=est_sex)) + 
    geom_density(size=2) + xlim(c(0,80)) +ggtitle(paste0(title," ", percentilex)) +
    scale_color_manual( values = c("male" = "light blue", "female" = "pink"), name="sex",
                        labels=c("male", "female")) 
  
  #(file.path(getwd(), "data-n-figs",Sys.Date()), showWarnings = FALSE)
  
  file.end <- ifelse(length(unique(post_data$est_sex))==2, "both-sexes.pdf","single-sex.pdf")
  ggsave(paste0(
    # getwd(),
    # "/data-n-figs/",
    # Sys.Date(),
    "/",
    percentilex,title,".",file.end), p, width = 12, height = 6, units="in")
  
  p<-NULL
  post_data <- NULL
  gc()
  
  
}
age_summary <- function(post_data, ...) {
  
  groupers<- quos(...)
  post_data %>%
    group_by(!!!groupers, est_sex) %>% 
    summarize(
      n=n(),
      sd=sd(est_age),
      mean=mean(est_age),
      mean.mrb = mean(mrb),
      mad=mad(est_age),
      median=median(est_age),
      IQR=IQR(est_age),
      ninety.perc.above = quantile(est_age,.1),
      nintey.perc.below = quantile(est_age,.9),
      prop.under.18.yrs = ecdf(est_age)(18)
    ) 
}



adult_prop_table <- function(post_data, estimated=T) {
  
  #either get ACTUAL or estimated ages and sex, depended on boolean "estimated"
  sex.expr <- ifelse(estimated==T, rlang::expr(est_sex), rlang::expr(sex))
  age.expr <- ifelse(estimated==T, rlang::expr(est_age), rlang::expr(age))
  
  post_data %>%
    group_by(!!sex.expr) %>%
    summarize(n.below.18= sum(!!age.expr<18), n.above.18= sum(!!age.expr>=18)) %>%
    column_to_rownames(as.character(sex.expr)) %>%
    prop.table() %>%
    rename_with(~gsub("n","p",.x)) #change n to p, for proportion
}



real_age_summary <- function(post_data, estimated=T, ...) {
  
  #either get ACTUAL or estimated ages and sex, depended on boolean "estimated"
  sex.expr <- ifelse(estimated==T, rlang::expr(est_sex), rlang::expr(sex))
  age.expr <- ifelse(estimated==T, rlang::expr(est_age), rlang::expr(age))
  
  
  groupers<- quos(...)
  post_data %>%
    group_by(!!!groupers, !!sex.expr) %>% 
    summarize(
      n=n(),
      sd=sd(!!age.expr),
      mean=mean(!!age.expr),
      mean.mrb = mean(mrb),
      mad=mad(!!age.expr),
      median=median(!!age.expr),
      IQR=IQR(!!age.expr),
      ninety.perc.above = quantile(!!age.expr,.1),
      nintey.perc.below = quantile(!!age.expr,.9),
      prop.under.18.yrs = ecdf(!!age.expr)(18)
    ) 
}

age_summary_function <- function(acol) {
  
  c(
    n=length(acol),
    sd=sd(acol),
    mean=mean(acol),
    mad=mad(acol),
    median=median(acol),
    IQR=IQR(acol),
    ninety.perc.above = quantile(acol,.1),
    nintey.perc.below = quantile(acol,.9),
    prop.under.18.yrs = ecdf(acol)(18)) %>% return()
  
}



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
  
#first, run a model with two sexes sharing an age distribution:
  cat("Model with 2 sexes, 1 age distribution: \n")
  model1age <-  db %>%
    run_model("2sex.1age", "2sex.1age.txt", pre.iters, pre.thin,
              n.chains=pre.n.chains,return.raw.samples=T, is.parallel = is.parallel)
  
#second, run a model with two sexes each with their own age distribution:
  cat("Model with 2 sexes, 2 age distributions: \n")
  model2age <- db %>%
    run_model("2sex.2age", "2sex.2age.txt", pre.iters, pre.thin,
              n.chains=pre.n.chains,return.raw.samples=T, is.parallel = is.parallel)

#thirdly, use these models to get your pseudopriors for the combined model
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
