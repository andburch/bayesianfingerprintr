
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesianfingerprintr

<!-- badges: start -->
<!-- badges: end -->

The goal of bayesianfingerprintr is to help archaeologists infer
demographic information about ancient potters who leave fingerprint
impressions on ceramic artifacts. It uses a Bayesian mixture modelling
approach, coupled with a data-driven understanding of how epidermal
ridge densities, sex, and age covary.

## Installation

This package requires the use of JAGS (Just Another Gibbs Sampler), and you'll need to install this software separately. 
Go to this page (https://sourceforge.net/projects/mcmc-jags/](https://sourceforge.net/projects/mcmc-jags/)) and hit "Download". The instructions should help you through the rest of the process, but if it doesn't work for you, you can go into the JAGS folder and find the specific version for your operating system.

Then you can install the development version of bayesianfingerprintr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("andburch/fingerprint-mixture-model")
```

## Example

This is a basic example, where we first make up some fake archaeological
mean ridge breadth (MRB) data. Each data point has an identification
number, a MRB measurement, and (we’re making this up here) maybe some
other information like where it was recovered from.

``` r
library(bayesianfingerprintr)

fake.data <- 
  data.frame(
  id.number = 1:100,
  mrb = c(rgamma(50,64,160), rgamma(50,2500,5000)),
  location = sample(c("palace", "domestic"), 100, replace = T)
  )
```

Then we want to make sure the MRB values are scaled to match the modern
reference data collected by Kralik and Novotny (2003):

``` r
fake.data$mrb <- fake.data$mrb*scaling_factor(fake.data$mrb)
```

Finally we can put these scaled MRB values into the model;

``` r
output <- 
  run_model(fake.data, #the data frame with the MRB values in a column named `mrb`
            "2sex.2age", #which variant of the model we want to run
            "2sex.2age.model.txt", #where to save the model specifications
            150, #number of MCMC iterations per chain (should be over 150,000 normally)
            5, #thinning factor
            F, F, 3, F, 100) #other options
```
