#kralik_data %>% filter(age<22) %>% {.$age*12 - .$mrb*614* 1.08108} %>% mean() %>% divide_by(112)

simulate_data <- function(Nsize = 100, prob_male = 0.8, MFmean_ages = c(20,15), MFsd_ages = c(10,10), errors=T, seed=1337) {
  
  library(qpcR, include.only = 'RMSE')
  
  fullfing <-read.csv(file=fingerprint_file_path) %>% as_tibble()
  
  kralik_data <-
    fullfing %>% 
    group_by(ID) %>%
    summarize(mrb=mean(RB), age=unique(age), sex=unique(sex)) 
  
  male_coef <-
    nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(kralik_data, sex==2)) %>% 
    coef()
  female_coef <-   
    nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(kralik_data, sex==1)) %>% 
    coef()
  
  set.seed(seed)
  
  prob_sex = c(prob_male, 1-prob_male)
  
  men_ages <- rgamma(Nsize,
                     (MFmean_ages[1]/MFsd_ages[1])^2,
                     MFmean_ages[1]/(MFsd_ages[1]^2))
  #hist(men_ages)
  
  
  women_ages <- rgamma(Nsize,
                       (MFmean_ages[2]/MFsd_ages[2])^2,
                       MFmean_ages[2]/(MFsd_ages[2]^2))
  #hist(women_ages)
  

  
  men_sd <- 
    nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(kralik_data, sex==2)) %>%
    qpcR::RMSE()
  men_errors <- rnorm(Nsize, 0, men_sd)/2
  
  women_sd <- 
    nls(mrb ~ SSasymp(age, Asym, r0, lrc), data=filter(kralik_data, sex==1)) %>%
    qpcR::RMSE()
  women_errors <- rnorm(Nsize, 0, women_sd)/2
  
  
  fake_dat <- 
    tibble(sex = sample(c("male","female"), Nsize, replace=TRUE, prob=prob_sex) %>% factor(levels = c("male","female"))) %>% 
    rowid_to_column() %>% 
    mutate(age = if_else(sex=="male",
                         men_ages[rowid],
                         women_ages[rowid])) %>% 
    mutate(mrb = case_when(
      sex == "male" ~ male_coef[1] + (male_coef[2] - male_coef[1])*exp(-1*exp(male_coef[3])*age),
      sex == "female" ~ female_coef[1] + (female_coef[2] - female_coef[1])*exp(-1*exp(female_coef[3])*age))
    ) 
  
  if(errors) fake_dat <- fake_dat %>% 
    mutate(mrb = if_else(sex=="male",
                         mrb + men_errors[rowid],
                         mrb + women_errors[rowid]))
  
  fake_dat <- fake_dat %>% 
    mutate(RD = 5*sqrt(2)/mrb,
           fake.RD = if_else(sex=="male",
                             floor(RD),
                             ceiling(RD))
           ) 
  
  return(fake_dat)
  

}

fowler_calculations <- function(db, male.cutoff=14.9, female.cutoff=14.9, shrinkage=0.8163196) {
  
  mrb_regression <- function(mrb, type=c("KAmod2","LC","LCmod", "KAmod","PM1", "random", "random2"), shrinkage) {
    case_when(
      type=="KAmod2" ~ (614*mrb-112*shrinkage)/12, #fowler's model
      type=="LC" ~ (1000*mrb-344.24)/11.94, # Loesch and Czyżewska (1972)
      type=="LCmod" ~ (mrb*1000*1.261 - 344.24)/11.94, #LC but not bad (by Kralik)
      type=="KAmod" ~ (614*1.08108*mrb-112)/12,  #Kralik's modified Kamp model
      type=="PM1" ~  52.18087*mrb- 7.89682,   #Kralik's data-generated model
      type=="random" ~ runif(1,5,60),  #randomly pick an age 5-60
      type=="random2" ~ rgamma(1,12.96,  0.72) #pick an age with mean 18, stdev 5
    ) %>% 
      return()
  }
  
  sex_estimation <- function(RD, male.cutoff, female.cutoff) { #these cutoffs are from Akiva, interestingly
    case_when(
      RD < male.cutoff ~ "male",
      RD >= male.cutoff & RD < female.cutoff ~ "indeterminate",
      RD >= female.cutoff ~ "female"
    ) %>% 
      # as.factor() %>% 
      return()
  }
  
  #add age estimates
  db <- bind_cols(db,
                  map_dfc(c("KAmod2","LC","LCmod", "KAmod","PM1", "random", "random2") %>% 
                            {set_names(.,paste0("est_age_",.))}, 
                          ~map_dbl(db$mrb, mrb_regression, type=., shrinkage=shrinkage)
                  )
  )
  db %>% 
    mutate(est_sex = map_chr(RD, sex_estimation,
                             male.cutoff=male.cutoff, female.cutoff=female.cutoff) %>% factor(levels=c("male","female")),
           est_sex_faked = map_chr(fake.RD, sex_estimation,
                                   male.cutoff=male.cutoff, female.cutoff=female.cutoff) %>% factor(levels=c("male","female")),
           est_sex_random = sample(c("male","female"),n(), replace=T) %>% factor(levels=c("male","female"))
           ) %>% 
    return()
  
  
  
}

RD_calc_verification <- function() {
  
  require("ggpubr")
  require("readxl")
  
  hama.data <- readxl::read_excel("C:/Users/aburchil/Downloads/andrew's spreadsheet.xlsx")
  
  hama.data <-
    hama.data %>% 
    mutate(mrb=breadth*10,
           calc.RD = 5*sqrt(2)/mrb,
           diff=`5x5 Ridges`-calc.RD) 
  
  print(sprintf("Correlation between measured RD and MRB: %.4f",
                cor(hama.data$`5x5 Ridges`, hama.data$mrb)))
  
  print(sprintf("Average difference between caclulated RD and measured RD: %.4f",
                mean(hama.data$diff)))
  
  print(sprintf("SD of difference between caclulated RD and measured RD: %.4f",
                sd(hama.data$diff)))
  
  print(sprintf("Mean *absolute* difference between caclulated RD and measured RD: %.4f",
                mean(abs(hama.data$diff))))
  
  
  p1<-
    ggplot(hama.data, aes(x=mrb, y=`5x5 Ridges`)) +
    labs(x="MRB", y="measured RD", title="MRB and RD are extremely highly correlated")+
    ylim(10,21.5)+
    geom_jitter(size=3, width=0, height=0.05) +
    geom_smooth(method="lm",formula=  y ~ I(1/x)) + theme_classic() + 
    ggpubr::stat_cor( aes(x=I(1/mrb), y=`5x5 Ridges`, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) 
  
  p2<-
    ggplot(hama.data, aes(y=calc.RD, x=`5x5 Ridges`)) +
    labs(y="RD calculated from MRB", x="measured RD", title="MRB explains 94% of variance in RD")+
    geom_jitter(size=3, height=0, width=0.05) +geom_smooth(method="lm") + theme_classic() + 
    ggpubr::stat_regline_equation()+ ylim(10,21.5)+
    ggpubr::stat_cor( label.x= 15,aes( label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) 
  
  ggpubr::ggarrange(p1,p2, ncol=2)
  
}

other_model_predictions <- function(db) {
  
  age_counts <- list (
    val.mean.age = ~mean(.),
    val.n.below.15= ~sum(.<15),
    val.n.15.to.20 = ~sum(.>=15 & . <=20),
    val.n.above.20= ~sum(.>20)
  )
  
  #get the total guess column estimates (they use their own sex estimates)
  actual <-
    db %>% select(sex,age) %>% 
    rename(est_sex ="sex") %>% 
    group_by(est_sex, .drop=F) %>% 
    summarize((across(where(is.numeric),age_counts)))
  
  randoms <-
    db %>% 
    select(contains("random")) %>% rename_at(vars(starts_with("est_age")), ~str_sub(.,9,-1)) %>% 
    rename(est_sex ="est_sex_random") %>% 
    group_by(est_sex, .drop=F) %>% 
    summarize((across(where(is.numeric),age_counts)))
  
  models <- 
    db %>% 
    select(-contains("random")) %>% 
    group_by(est_sex_faked) %>% select(starts_with("est_age_")) %>% 
    rename_at(vars(starts_with("est_age")), ~str_sub(.,9,-1)) %>% 
    rename(est_sex ="est_sex_faked") %>% 
    group_by(est_sex, .drop=F) %>% 
    summarize((across(where(is.numeric),age_counts))) %>% 
    full_join(actual,., by="est_sex") %>% 
    full_join(randoms, by="est_sex")
  
  model_tibble <-
    models %>% 
    pivot_longer(cols = contains("_val."),values_to="n" ) %>%
    separate(name, c("model","prediction"), sep="_") %>%
    unite("predictions", c(est_sex,prediction)) %>%
    {split(.,.$model)} %>% 
    #turn them into percentages and rename everything
    map(.,
        ~{transmute(., 
                    predictions = str_replace(predictions, "val.n.","val.perc."),
                    predictions = str_replace(predictions, "_val.","_"),
                    value = ifelse(str_detect(predictions, "mean.age"),
                                   n, n/sum(n[str_which(predictions, "_perc.")]))
        ) %>% 
            arrange(predictions)}
    ) %>%
    #add a few more things
    map(.,
        ~{perc.female <- filter(., str_detect(predictions, "^female_perc"))$value %>% sum();
        mean.age <- 
          weighted.mean(c(filter(., str_detect(predictions, "mean.age"))$value),
                        c(perc.female, 1-perc.female), na.rm=T);
        add_row(., predictions = "total_mean.age", value=mean.age) %>% 
          add_row(., predictions = "perc_female", value=perc.female)}) %>% 
 
    t() %>% as_tibble() %>% rename_with(~paste0("model_",.))
  
  return(model_tibble)
  
  
}

run_single_model_ALOT <- function(db,iters=1500, thin=3) {
  age_counts <- list (
    val.mean.age = ~mean(.),
    val.n.below.15= ~sum(.<15),
    val.n.15.to.20 = ~sum(.>=15 & . <=20),
    val.n.above.20= ~sum(.>20)
  )
  
  statz <- 
    db %>% select(mrb) %>% 
    run_single_model(iters=iters, thin=thin) %>% 
    select(est_sex, bayes=est_age) %>% 
    group_by(est_sex, .drop=F) %>% 
    summarize((across(where(is.numeric),age_counts))) %>% 
    pivot_longer(cols = contains("_val."),values_to="n" ) %>%
    separate(name, c("model","prediction"), sep="_") %>%
    unite("predictions", c(est_sex,prediction)) %>%
    {split(.,.$model)} %>% 
    #turn them into percentages and rename everything
    map(.,
        ~{transmute(., 
                    predictions = str_replace(predictions, "val.n.","val.perc."),
                    predictions = str_replace(predictions, "_val.","_"),
                    value = ifelse(str_detect(predictions, "mean.age"),
                                   n, n/sum(n[str_which(predictions, "_perc.")]))
        ) %>% 
            arrange(predictions)}
    ) %>%
    #add a few more things
    map(.,
        ~{perc.female <- filter(., str_detect(predictions, "^female_perc"))$value %>% sum();
        mean.age <- 
          weighted.mean(c(filter(., str_detect(predictions, "mean.age"))$value),
                        c(perc.female, 1-perc.female), na.rm=T);
        add_row(., predictions = "total_mean.age", value=mean.age)}) %>% 
    t() %>% as_tibble() %>% rename_with(~paste0("model_",.)) 
  
  db<-NULL
  
  
  return(statz)
  
  
}

create_many_scenarios <- function() {
  #this IS the HPC  model that I actually ran
  #lol forgot to use multi-cores
  set.seed(111)
  
  huge.db <-  expand.grid(Nsize = c(50,100,200), 
                          prob_male = seq(.1,.9,.2),
                          meanM= c(15,20,30),
                          meanF = c(15,20,30),
                          sdM = c(2,8),
                          sdF= c(2,8)) %>% 
    as_tibble() %>% 
    mutate(seed=sample.int(1000,n()))
  
  
  huge.db <-
    huge.db %>% rowwise() %>% 
    mutate(simulated_data = list(simulate_data(Nsize, prob_male, c(meanM, meanF), c(sdM, sdF), seed))) %>% ungroup()
  
  
  huge.db <-
    huge.db %>% rowwise() %>% 
    mutate(other_model_estimates = list(fowler_calculations(simulated_data,14.9,14.9))) %>% ungroup()
  
  huge.db<- 
    huge.db %>% rowwise() %>%
    mutate(other_model_predictions(other_model_estimates)) %>% ungroup()
  
  # huge.db <-
  #   huge.db %>% rowwise() %>%
  #   mutate(run_single_model_ALOT(simulated_data, 500000,5)) %>% ungroup()
  # 
  holder<-NULL
  for (i in 1:nrow(huge.db)) {
    print(huge.db[i,1])
    holder <- rbind(holder,run_single_model_ALOT(huge.db[i,]$simulated_data[[1]],5000,5))
    
  }
  huge.db<-cbind(huge.db,holder) %>% as_tibble()
  
  
  saveRDS(huge.db,"huuugedb.RDS")
  
  
}


#heatmap of model success
plot_models_success <- function(dbz,sex.est=T) {
  
  #this adds age-based metrics irrespective of sex
  add_sexfree_age_metrics <- function(db) {
    
    summo <- function(db.subset) {
      
      db.subset %>% 
        summarize_at(vars(values, starts_with("model_")), sum) %>% #sum up the male and female percents
        cbind(db.subset %>% select(!starts_with("model_"), -values),.) %>% #add them back to the normal data
        .[1,] %>% #select only one row
        mutate(measures=paste0("total", str_match(measures,"_perc.+"))) #rename it
      
    }
    db<- db %>% filter(str_detect(measures, "perc.15.")) %>% summo() %>% rbind(db,.)
    db<- db %>% filter(str_detect(measures, "perc.above.20")) %>% summo() %>% rbind(db,.)
    db<- db %>% filter(str_detect(measures, "perc.below.15")) %>% summo() %>% rbind(db,.)
    return(db)
    
  }
  
  #print/make heatmap
  heatmappy <- function(db) {
    
    db %>% 
      group_by(measures) %>% 
      #remove extra padding etc around model names
      select(starts_with("model_")) %>% rename_at(vars(starts_with("model_")), ~str_sub(.x,7,-8)) %>% 
      summarize_all(mean, na.rm=T) %>% ungroup() %>%
      #add "measures" column to matrix rownames and remove from matrix itself
      {set_rownames(as.matrix(.[,-1]), .[[1]])} %>% 
      #reorder columns based on means
      {.[,order(colMeans(., na.rm=TRUE))]} -> matvals
    
    order.values<- matvals %>% 
      #scale the rows so that 
      {t(scale(t(.)))} %>% 
      #reorder columns based on means
      {.[,order(colMeans(., na.rm=TRUE))]} %>% 
      colnames()
      
    color.values <- matvals %>% 
      #scale the rows so that it looks nice
      {t(scale(t(.)))} %>% 
      as.data.frame() %>% 
      rownames_to_column(var="measures") %>% 
      pivot_longer(-c(measures), names_to="models", values_to = "colors") %>%
      mutate(models=factor(models, levels=order.values))
    
    actual.values <-matvals %>% 
      as.data.frame() %>% 
      rownames_to_column(var="measures") %>% 
      pivot_longer(-c(measures), names_to="models", values_to = "actuals") %>%
      mutate(actuals=if_else(str_detect(measures,"perc"), actuals*100,actuals),
             models=factor(models, levels=order.values))
    
    full_join(actual.values,color.values) %>%
      mutate(measures= measures %>%
               str_replace("(.+)(perc\\.)","% \\1") %>%
               str_replace("ale","ales") %>% 
               str_replace_all("[_\\.]"," ") %>% 
               str_replace("perc","%")
      ) ->x
    xx<<-x
      ggplot(xx,aes(x=models,y=measures,fill=actuals, label=round(actuals, 2))) +
      geom_tile() +
      scale_fill_gradient2(low="dark green",high="red", midpoint=15) +
      geom_text()+theme_classic() + guides(fill = "none") %>% 
      return()
      
   
      #reorder ROWS based on means
      #{.[order(rowMeans(., na.rm=TRUE)),]} %>% 
      
      # heatmap(Colv = NA, Rowv = NA, scale="row", margins=c(6,8), 
      #         main="Difference between model estimates and actual measures",
      #         col=colorRampPalette(colors = c("green", "white", "red"))(100)) %>% 
      # print()
    
    
    
  }
  
  dbz <- dbz %>%
    select(-simulated_data,-other_model_estimates) %>% 
    unnest(cols=starts_with("model"), names_sep = "..") %>% 
    group_by(model_age..predictions) %>%  select(-ends_with("predictions")) %>% ungroup() %>%
    rename("measures"=model_age..predictions,"values"=model_age..value) %>% 
    #filter(measures=="perc_female") %>% 
    filter(!is.nan(values)) %>%
    group_by(seed) %>% group_split() %>% map_df(add_sexfree_age_metrics) %>%  #add age metrics that don't have sex
    mutate_at(vars(starts_with("model_")), list(~abs(values-.))) #find differences

  #
  if(sex.est==T) {
    dbz %>%  filter(!str_detect(measures, "total_perc.+")) %>% heatmappy() %>% return()
  } else {
    dbz %>%  filter(str_detect(measures, "total_.+")) %>% heatmappy()%>% return()
  }
  
  
}

#====Get sex threshold=======


#based on Kraliks data, this is what the threshold should be:
sex.threshold.model <- glm(sex ~ RD, 
                           data=kralik_data %>% mutate( RD = 5*sqrt(2)/mrb, sex=sex-1) %>% 
                             filter(age>18), family="binomial") 

sex.threshold <- -coef(sex.threshold.model)[1] / coef(sex.threshold.model)[2]


#====Get big db with ran-scenarios========

add_perc_female <-function(db) {
  bind_rows(db,
            list(predictions="perc_female",
                 value=subset(db,str_detect(predictions,"^female_perc.+"))$value %>% sum()))
  
}
db1<-readRDS("C:/Users/aburchil/Downloads/huuugedb (1).RDS")
db2<-readRDS("C:/Users/aburchil/Downloads/huuugedb2.RDS")
db3<-readRDS("C:/Users/aburchil/Downloads/unfixednow.RDS") #old gamma, unfixed

db4<-readRDS("C:/Users/aburchil/Downloads/unfixednow (1).RDS") #new gamma, unfixed

db5<-readRDS("C:/Users/aburchil/Downloads/unfixednow (2).RDS") #new gamma, fixed


db1$model_bayes1age <- db2$model_bayes
db1$model_unfixed <-db3$model_bayes
db1$model_unfixedgamma <-db4$model_bayes
db1$model_fixedgamma <-db5$model_bayes



db1 <-
  db1 %>%
  mutate_at(vars(starts_with("model_")), ~map(.,add_perc_female))

plot_models_success(
  db1 %>% select(-model_LCmod,-model_LC,-model_KAmod,-model_PM1) %>%
    rename(model_random="model_bozo", model_random2="model_bozo2")
  )


  
  
  
###### 
  # 
  # db1 %>%   unnest(cols=starts_with("model"), names_sep = "..") %>% 
  #   group_by(model_age..predictions) %>%  select(-ends_with("predictions")) %>% ungroup() %>%
  #   rename("measures"=model_age..predictions,"values"=model_age..value) %>% 
  #   select(-simulated_data,-other_model_estimates) %>% 
  #   mutate_at(vars(starts_with("model_")), list(~abs(values-.))) %>% 
  #   pivot_wider(names_from = measures, values_from = model_bayes..value, id_cols=1:8 ) %>% 
  #   
  #   lm(female_mean.age ~ Nsize*prob_male*meanM*meanF*sdM*sdF, data=.) %>% summary()
  # 
  # 
  
  

  
  
  #################
  
  kralik_data %<>% rowwise() %>% mutate(RD = 5*sqrt(2)/mrb,
                          fake.RD = if_else(sex=="male",
                                            floor(RD),
                                            ceiling(RD)),
                          sex=ifelse(sex==2,"male","female"),
                          sex=as.factor(sex)
  )
  
  hey<-
  mutate(tibble(type="kralik"),
         simulated_data = list(kralik_data )) %>% 
    rowwise() %>% 
         mutate(other_model_estimates = list(fowler_calculations(kralik_data,14.9,14.9)),
         other_model_predictions(other_model_estimates)) %>% 
    rowwise() %>% 
    mutate(run_single_model_ALOT(simulated_data, 150000,5))
  
  hey$model_bayes[[1]] %>% add_perc_female() -> hey$model_bayes[[1]] 
  hey$seed<-1
  
  plot_models_success(hey)
  
  
  
  hey %>%   unnest(cols=starts_with("model"), names_sep = "..") %>% 
    group_by(model_age..predictions) %>%  select(-ends_with("predictions")) %>% ungroup() %>%
    rename("measures"=model_age..predictions,"values"=model_age..value) %>% 
    filter(measures=="percent") %>% select(-simulated_data,-other_model_estimates) %>% 
    mutate_at(vars(starts_with("model_")), list(~abs(values-.))) %>% select(starts_with("model")) %>%
    colMeans(na.rm=T) %>% sort() %>% barplot()
  






  


