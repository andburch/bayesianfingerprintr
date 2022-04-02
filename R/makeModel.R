#' Create the JAGS model file
#'
#' This writes the desired fingerprint mixture model into a separate file to be read into JAGS. It also returns the nodes that need to be monitored during the MCMC run.
#'
#' @param model.opt a string representing which version of the model you want to run
#' @param filename the name of the text document that the model will be written to
#' @param pseudopriors when running a "super model," this is where you pass the posterior distributions of the two models to then be used as pseudopriors
#'
#' @return a list of the nodes (parameters) that should be monitored during the JAGS MCMC run
#'
#' @author Andrew Burchill, \email{andrew.burchill@asu.edu}
#'
#' @examples
#' print("hello")
#'
#' @import dplyr
#' @rawNamespace import(stats, except = c(lag, filter))


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


ref_data <-
  ref_data %>%
  as_tibble() %>%
  group_by(.data$ID) %>%
  summarize(mrb=mean(.data$RB), age=unique(.data$age), sex=unique(.data$sex))

male_coef <-
  nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=dplyr::filter(ref_data, .data$sex==2)) %>%
  coef()
female_coef <-
  nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=dplyr::filter(ref_data, .data$sex==1)) %>%
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
