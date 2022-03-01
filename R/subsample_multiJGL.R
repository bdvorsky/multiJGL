

#' @title Function for running multiJGL algorithm across subsamples and different sparsity levels
#'
#' @param node.covariates An nxp dimensional matrix of p covariates measured over n samples.
#' @param grouping.factor A grouping factor for creating observational classes.
#'
#' @param lambda.seq The range of different lambda values
#' @param by_value Specifies the density of lambda grid
#' @param num_repetitions The number of subsample analyses
#' @param penalty.lin Specify "fused" or "group" penalty type for the linear JGL algorithm.
#' @param penalty.nonlin Specify "fused" or "group" penalty type for the nonlinear JGL algorithm.
#' @param lin_lambda1 The l1-penalty parameter for the linear JGL to regulate within group network densities
#' @param lin_lambda2 The l1-penalty parameter for the nonlinear JGL.
#'
#' @param nonlin_lambda1 The fusion penalty parameter for the linear JGL.
#' @param nonlin_lambda2 The fusion penalty parameter for the nonlinear JGL.
#' @param tol.linear Convergence criterion for the linear part (see the JGL package for details).
#' @param tol.nonlinear Convergence criterion for the nonlinear part.
#' @param subsample.pseudo_obs Should the subsampling procedure be used over the pseudo-observations.
#' @param omit.rate An integer: Omit rate for the subsampling pcocedure between 2L and 5L

#'
#'
#' @examples print("")
#' @export

subsample_multiJGL<- function(node.covariates = node.covariates,
                              grouping.factor = grouping.factor,

                              lambda.seq,
                              by_value,
                              num_repetitions = 10,
                              penalty.lin = "fused",
                              penalty.nonlin = "fused",
                              lin_lambda1 = lambda.seq1[j],
                              lin_lambda2 = 0.025,
                              nonlin_lambda1 = lambda.seq1[j],
                              nonlin_lambda2 = 0.025,
                              tol.linear = 1e-05,
                              tol.nonlinear = 1e-05,
                              subsample.pseudo_obs = FALSE,
                              omit.rate = 2L){

  #helper-function for creating and naming list objects
  create.list <- function(num.of.objects, object.names){
    list <- vector(mode = "list", length = num.of.objects)
    names(list) <- object.names
    return(list)
  }


  num.of.classes <- length(unique(grouping.factor))
  obs.classes <- vector(mode = "list", length = num.of.classes)

  for(i in 1:num.of.classes){
    obs.classes[[i]] <- node.covariates[which(as.numeric(as.factor(grouping.factor)) == i),]
  }

  X <- obs.classes

  #Set up a grid of different lambda1 (JGL algorithm, fused penalty) values.
  lambda.seq1 <- lambda.seq
  #The number of subsamples and networks
  num_repetitions <- num_repetitions;
  num_networks <- length(X)

  #Make a lambda and network specific list of estimated networks over subsamples
  lin_lambda_repetitions_list <- nonlin_lambda_repetitions_list <-
       lapply(create.list(length(lambda.seq1),
                          paste("lambda=", lambda.seq1, sep = "")),
  #Each object in layer (2) represents results from different subsample analyses.
        function(x){list <- create.list(num_repetitions,
                            paste("subsample", seq(1:num_repetitions), sep = "")) %>%
  #Each object within layer (3) contains subsample-specific linear networks.
          lapply(function(x)  create.list(num_networks,
                              paste("Network", seq(1:num_networks), sep = "")))
          return(list)})

  for(j in 1:length(lambda.seq1)){
    for(l in 1:num_repetitions){
      subsample_X <- vector(mode = "list", length = length(X))
      for(k in 1:length(X)){
        subsample_X[[k]] <- X[[k]][-sample(nrow(X[[k]]),ceiling(0.2*nrow(X[[k]])), replace = FALSE),]
        nw <-
          multiJGL(node.covariates = subsample_X,
                         grouping.factor = grouping.factor,
                         penalty.lin = "fused",
                         penalty.nonlin = "fused",
                         lin_lambda1 = lambda.seq1[j], lin_lambda2 = 0.025,
                         nonlin_lambda1 = lambda.seq1[j],nonlin_lambda2 = 0.025,
                         tol.linear = 1e-05,
                         tol.nonlinear = 1e-05,
                         subsample.pseudo_obs = FALSE,
                         omit.rate = 2L)

    lin_lambda_repetitions_list[[j]][[l]] <- nw[[1]]$theta
              names(lin_lambda_repetitions_list[[j]][[l]]) <-
              paste("Network", seq(1:length(X)), sep = "")


    nonlin_lambda_repetitions_list[[j]][[l]] <- nw[[2]]$theta
              names(nonlin_lambda_repetitions_list[[j]][[l]]) <-
              paste("Network", seq(1:length(X)), sep = "")
      }
    }
  }
}
