
#' Title: A dynamic network model estimator from continuous grouping variables
#' @param response A continuous response vector.
#' @param node.covariates An nxp dimensional matrix of p covariates measured over n samples.
#' @param rate Defines the step length for the segment slidings.
#' @param segment.width An empirical distribution segment width.
#' @param penalty.lin Specify "fused" or "group" penalty type for the linear JGL algorithm.
#' @param penalty.nonlin Specify "fused" or "group" penalty type for the nonlinear JGL algorithm.
#' #See the explanation for the fused and group penalties in the JGL package
#' #The original JGL CRAN repository: https://CRAN.R-project.org/package=JGL
#' #The following penalty parameters are given in pairs to --
#' #separately assign the amount of regularizations for linear and nonlinear parts
#' @param lin_lambda1 The l1-penalty parameter for the linear JGL to regulate within group network densities
#' @param lin_lambda2 The l1-penalty parameter for the nonlinear JGL.
#'
#' @param nonlin_lambda1 The fusion penalty parameter for the linear JGL.
#' @param nonlin_lambda2 The fusion penalty parameter for the nonlinear JGL.
#' @param tol.linear Convergence criterion for the linear part (see the JGL package for details).
#' @param tol.nonlinear Convergence criterion for the nonlinear part.
#' #Subsampling procedure over pseudo-observations if the number of observations is large already in the original sets.
#' @param subsample.pseudo_obs Should the subsampling procedure be used over the pseudo-observations.
#' @param omit.rate An integer: Omit rate for the subsampling pcocedure between 2L and 5L

#' @import whitening
#' @import JGL
#' @importFrom CVglasso CVglasso
#' @importFrom matrixStats rowDiffs
#' @importFrom stats cov
#' @importFrom utils combn
#' @importFrom stats quantile
#' @examples  print("net <- dynamic_NLJGL(node.covariates, grouping.factor)")
#' @export
dynamicMultiJGL<- function(node.covariates = node.covariates, response = response,
                           segment.width = 1/2,
                           rate = 1/10,
                           penalty.lin = "fused",
                           penalty.nonlin = "fused",
                           lin_lambda1 = 0.025, lin_lambda2 = 0.01,
                           nonlin_lambda1 = 0.025,nonlin_lambda2 = 0.01,
                           tol.linear = 1e-05,
                           tol.nonlinear = 1e-05,
                           subsample.pseudo_obs = FALSE,
                           omit.rate = 2L){

  #Generate the rate and width of the moving segment
  low_seq <- seq(0, segment.width, by = rate)
  high_seq <- seq(segment.width, 1.0, by = rate)

  #This algorithm is built upon the linear JGL R ´
  #CRAN repository: https://CRAN.R-project.org/package=JGL

  #Check if the input dataset contains missing values
  if(anyNA(node.covariates)) stop("NAs in node.covariates data not allowed")


  #Set up the parameters
  numb.objects <- length(low_seq)
  obs.classes <- vector(mode = "list", length = numb.objects)
  whitened_data_matrices <- vector(mode = "list", length = numb.objects)
  pseudo_obs <- vector(mode = "list", length = numb.objects)

  #The whitening procedure to remove linear relationships
  for(i in 1:numb.objects){
    obs.classes[[i]] <- node.covariates[which(response <= stats::quantile(response, high_seq[i]) & response >= quantile(response, low_seq[i])),]

    Ctest <- CVglasso::CVglasso(X = obs.classes[[i]], S = cov(obs.classes[[i]]),
                      #Apply a small penalty parameter lambda value just to get a non-singular input for the ZCA-cor method.
                      nlam = 10, lam.min.ratio = 0.0,
                      lam = 0.001, diagonal = FALSE, path = FALSE,
                      tol = 1e-04,crit.cv = "AIC",
                      maxit = 10000)

    S <- Ctest$Sigma
    #Calculate the whitening matrix based on the covariance matrix
    W.ZCAcor = whiteningMatrix(S, method="ZCA-cor")
    #Apply the ZCA-cor whitening procedure
    whitened_data_matrices[[i]] = tcrossprod(as.matrix(obs.classes[[i]]), W.ZCAcor)
    whitened_domain <- scale(whitened_data_matrices[[i]])
    whitened_domain <- whitened_domain[order(whitened_domain[,1]),]
    tr_whitened_domain <- t(whitened_domain)

    #Pseudo-observation generating step
    #Enumerate all pairwise combinations between sample indexes
    cols <- combn( ncol(tr_whitened_domain), 2)
    tripleSums <- apply( cols, 2, function(z) rowDiffs(tr_whitened_domain[,z]))
    whitened_domain <- t(tripleSums)
    #Generate kernel specific pseudo-observation sets
    pseudo_obs[[i]] <- scale(sqrt(abs(whitened_domain)))

    colnames(pseudo_obs[[i]]) <- colnames(node.covariates)
  }

  #Estimate the linear and nonlinear network structures with the original JGL-algorithm

  lin_fgl.results = JGL(Y=obs.classes,lambda1=lin_lambda1, lambda2=lin_lambda2,
                        penalty = penalty.lin, return.whole.theta = TRUE,
                        tol = tol.linear)

  #This part subsamples the pseudo-observation set
  if(subsample.pseudo_obs == TRUE){

    if(!(is.integer(omit.rate))) stop("omit.rate must be an integer (2L, 3L, 4L or 5L)")

    if(omit.rate < 2 | omit.rate > 5) stop("omit.rate must be between 2 and 5")
    sub.pseudosample <-  lapply(pseudo_obs, function(class){

      return(class[c(TRUE,rep(FALSE, omit.rate-1)),])})

    nonlin_fgl.results = JGL(Y=sub.pseudosample, penalty=penalty.nonlin,
                             lambda1=nonlin_lambda1,
                             lambda2=nonlin_lambda2,
                             return.whole.theta = TRUE,
                             tol = tol.nonlinear)
  } else {
    nonlin_fgl.results = JGL(Y=pseudo_obs,lambda1=nonlin_lambda1,lambda2=nonlin_lambda2,
                             penalty = penalty.nonlin, return.whole.theta = TRUE,
                             tol = tol.nonlinear)
  }


  results <- list(linear_networks = lin_fgl.results, nonlinear_networks = nonlin_fgl.results)

  return(results)
}
