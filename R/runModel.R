#' Run the MCMC sampler for the Bayesian mixture model
#'
#' This function, as with all Bayesian models, will probably take a long time to run.
#'
#' @param db the data frame of observed MRB values (contained in a column named 'mrb')
#' @param model.opt select which version of the model to run
#' @param filename the name of the file that the JAGS model will be saved to
#' @param iters the number of iterations of the MCMC sampler
#' @param thin thinning rate in the MCMC sampler. Must be a positive integer.
#' @param is.parallel if parallel=TRUE, specify the number of CPU cores used. Defaults to total available cores or the number of chains, whichever is smaller.
#' @param DIC option to report the deviance information criterion and the estimated number of parameters (pD). (perhaps useful for model selection)
#' @param n.chains the number of chains in the MCMC sampler
#' @param return.raw.samples whether the raw MCMC posterior samples be returned or a tidy database of samples linked to the actual input data
#' @param n.adapt number of iterations to run in the JAGS adaptive phase. The default is NULL, which will result in the function running groups of 100 adaptation iterations (to a max of 10,000) until JAGS reports adaptation is sufficient. If you set n.adapt manually, 1000 is the recommended minimum value
#' @param n.burnin number of iterations at the beginning of the chain to discard (i.e., the burn-in). Does not include the adaptive phase iterations.
#' @param fix.biggest if TRUE, the data point with the largest MRB will be fixed as male. This can help if there are issues of model identifiability.
#' @param pseudopriors (optional) this is for passing pseudopriors into the model when running a "super model"
#'
#' @return depending on the `return.raw.samples` parameter, either a MCMC list object or a data frame with the posterior distribution joined with the original data set.
#' @export
#'
#' @author Andrew Burchill, \email{andrew.burchill@asu.edu}
#'
#' @examples
#' \dontrun{
#' db = data.frame(mrb=c(rgamma(50,64,160),
#'                       rgamma(50,2500,5000)),
#'                location=sample(c("palace", "domestic"), 100, replace=T)
#'                )
#' output = run_model(db, "2sex.2age", 15000, 5, F, F, 3, F)
#'
#' }
#'
#' @import stats dplyr

run_model <- function(
  db, model.opt = c("2sex.2age","1sex.1age","2sex.1age","1sex.2age", "compare1v2age", "compare1v2ageNEXT"), filename="bayesmodel.txt",
  iters, thin, is.parallel=F, DIC=F, n.chains=10, return.raw.samples=F,
  n.adapt=3000, n.burnin=100, fix.biggest=F,pseudopriors=NULL){

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
