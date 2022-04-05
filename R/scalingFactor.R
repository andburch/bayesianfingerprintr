#' Align MRB assemblages to the modern reference dataset
#'
#' Interassemblage variation makes comparing mean ridge breadth (MRB) measurements
#' difficult. This function scales a vector of MRB values to the modern
#' Kralik and Novotny (2003) reference data set.
#'
#'
#' @param mrbs the vector or column of your collected MRB measurements
#' @param x what percentile should serve as the alignment point? (defaults to 95th percentile)
#' @param type are you aligning your data to a subset of the reference data? Options are currently \code{"all","males","females","adult males","adult females","adults"}. We highly suggest sticking with the default unless there is strong evidence otherwise.
#' @param adult_age what age do you consider to be the threshold for adulthood? (only used when scaling to a subset of the reference data)
#'
#' @return the scaling factor that your data should be multiplied by in order to be analagous with modern MRB data
#' @export
#'
#' @examples
#' fake.data <-
#' data.frame(
#'  id.number = 1:100,
#'   mrb = c(rgamma(50,64,160), rgamma(50,2500,5000)),
#'   location = sample(c("palace", "domestic"), 100, replace = TRUE)
#' )
#'
#' fake.data$mrb <- fake.data$mrb*scaling_factor(fake.data$mrb)
#'
#' @import dplyr

scaling_factor <- function(mrbs, x=0.95, type="all", adult_age=15) {
  if (!(type %in% c("all","males","adult males","adult females","adults", "females"))) stop("You need to select an appropriate subsetting option.")
  db = switch(type,
              "all" = ref_data,
              "adults" = dplyr::filter(ref_data, .data$age >= adult_age),
              "males" = dplyr::filter(ref_data, .data$sex == 2),
              "females" = dplyr::filter(ref_data, .data$sex == 1),
              "adult males" = dplyr::filter(ref_data, .data$age >= adult_age, .data$sex == 2),
              "adult females" = dplyr::filter(ref_data, .data$age >= adult_age, .data$sex == 1)
              )
  quantile(db$mrb,x)/quantile(mrbs,x)
}
