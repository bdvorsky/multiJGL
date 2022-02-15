#' Anomaly checking for network node variables
#'
#' @param covariates A numeric sample matrix with n rows (all samples) and p columns of network covariates
#' @param groups A character vector indicating a group in which each sample belongs to.
#'
#' @return The first table provides information about the empirical distributions of each variable:
#'  Min, max, skewness, and kurtosis. The second table shows group-specific sample sizes.
#'
#' @examples print("anomaly_check(covariates)")
#' @importFrom e1071 kurtosis
#' @importFrom e1071 skewness
#' @export

data_check <- function(covariates, groups){

  message("Information about the variable-specific empiricical distributions:")

 return(data.frame(min = apply(covariates, 2, function(x) min(x)),
             max = apply(covariates, 2, function(x) max(x)),
             skewness = apply(covariates, 2, function(x) skewness(x)),
             exc_kurtosis = apply(covariates, 2, function(x) kurtosis(x))) %>%
          round(., 3) %>%
    mutate("warning_I" = case_when(
      skewness < -1 ~ "skwns < -1",
      skewness >  1 ~ "skwns > 1",
      TRUE ~ "-")) %>%
      mutate("warning_II" = case_when(
        exc_kurtosis <  -1 ~ "krts < -1",
        exc_kurtosis >  5 ~ "krts > 5",
        TRUE ~ "-")))

  return(table(groups))


}
