#' Title dynamic plot creation function
#'
#' @param networks Estimated dynamic network structure
#' @param display.networks Whether LINEAR, NONLINEAR or COMBINED networks are considered
#' @param display.unique Whether all structures or only group specific unique structures are plotted
#' @param main Title for the animation
#' @param xlab xlab for the animation
#' @param output.mode output mode
#' @param ... additional arguments for the render.d3movie function used to create animations
#' @importFrom network as.network
#' @importFrom networkDynamic networkDynamic
#' @importFrom ndtv render.d3movie
#' @return Dynamic output, e.g., as htmlWidget
#' @examples print("dynamicMultiJGLplot(dynamicMultiJGL(x))")
#' @export
dynamicMultiJGLplot <- function(networks,
                                display.networks = c("LINEAR"),
                                display.unique = "FALSE",

                                main = NULL,
                                xlab=NULL,

                                output.mode = "htmlWidget", ...){

  if(is.null(main)){
    cat("Argument main is null:  Default main is used")
    main <- "Observational class specific network structures"
  }
  if(is.null(xlab)){
    cat("Argument xlab is null: Default xlab is used\n")
    xlab <- "Observational classes"
  }


  #Initialize objects for the algorithm
  desired_length <- length(networks[[1]]$theta)
  lin_fgl.results <- networks[[1]]
  nonlin_fgl.results <- networks[[2]]


  linear_network_list <- vector(mode = "list", length = desired_length)
  nonlinear_network_list <- vector(mode = "list", length = desired_length)
  combined_network_list <- vector(mode = "list", length = desired_length)
  combined_mat_list <- vector(mode = "list", length = desired_length)



  for(i in 1:desired_length){

    linear_network_list[[i]] <- as.network(lin_fgl.results$theta[[i]])
    nonlinear_network_list[[i]] <- as.network(nonlin_fgl.results$theta[[i]])
    combined_network_list[[i]] <- as.network(abs(nonlin_fgl.results$theta[[i]]) +
                                               abs(lin_fgl.results$theta[[i]]))
    combined_mat_list[[i]] <- abs(nonlin_fgl.results$theta[[i]]) +
      abs(lin_fgl.results$theta[[i]])

  }

  # convert a list of networks into networkDynamic object

  if(display.unique == "FALSE"){

    if(display.networks == "LINEAR"){
      tnet<-networkDynamic(network.list=linear_network_list)
      print(render.d3movie(tnet,
                           main = main,
                           xlab=xlab, displaylabels =TRUE ,
                           output.mode = output.mode, ...))
    }
    if(display.networks == "NONLINEAR"){
      nonlintnet<-networkDynamic(network.list=nonlinear_network_list)
      print(render.d3movie(nonlintnet,
                           main = main,
                           xlab=xlab, displaylabels =TRUE,
                           output.mode = output.mode, ...))
    }
    if(display.networks == "COMBINED"){
      combinedtnet<-networkDynamic(network.list=combined_network_list)
      print(render.d3movie(combinedtnet,
                           main = main,
                           xlab=xlab, displaylabels =TRUE ,
                           output.mode = output.mode, ...))
    }
  }

  #Extract the unique network structures
  if(display.unique == "TRUE"){

    unique_struct <- function(adj_mat_list, i){

      #Initialize the complementary adjacency matrix
      complementary_list <- adj_mat_list
      reference_adj_mat <- complementary_list[[i]]
      complementary_list[[i]] <- matrix(0, ncol(adj_mat_list[[1]]),
                                        ncol(adj_mat_list[[1]]))
      #Compute the element wise sum from the above list of matrices
      complementary_adj_mat <- Reduce('+', complementary_list)
      complementary_adj_mat[which(complementary_adj_mat != 0)] <- 1
      reference_adj_mat[which(reference_adj_mat != 0)] <- 1
      #Compute the unique structures
      diff <- reference_adj_mat - complementary_adj_mat

      diff[which(diff == -1)] <- 0
      colnames(diff) <- colnames(adj_mat_list[[1]])
      return(diff)

    }


    lin_unique_structure_network_list <- vector(mode = "list",
                                                length = desired_length)
    nonlin_unique_structure_network_list <- vector(mode = "list",
                                                   length = desired_length)
    combined_unique_structure_network_list <- vector(mode = "list",
                                                     length = desired_length)



    for(i in 1:desired_length){

      lin_unique_structure_network_list[[i]] <- as.network(unique_struct(lin_fgl.results$theta,i))
      nonlin_unique_structure_network_list[[i]] <- as.network(unique_struct(nonlin_fgl.results$theta,i))
      combined_unique_structure_network_list[[i]]<- as.network(unique_struct(combined_mat_list,i))
    }


    if(display.networks == "LINEAR"){
      tnet<-networkDynamic(network.list=lin_unique_structure_network_list)
      print(render.d3movie(tnet,
                           main = main,
                           xlab= xlab, displaylabels =TRUE ,
                           output.mode = output.mode, ...))
    }
    if(display.networks == "NONLINEAR"){
      nonlintnet<-networkDynamic(network.list=nonlin_unique_structure_network_list)
      print(render.d3movie(nonlintnet,
                           main = main,
                           xlab=xlab ,displaylabels =TRUE ,
                           output.mode = output.mode, ...))
    }
    if(display.networks == "COMBINED"){
      combinedtnet<-networkDynamic(network.list=combined_unique_structure_network_list)
      print(render.d3movie(combinedtnet,
                           main = main,
                           xlab= xlab,displaylabels =TRUE ,
                           output.mode = output.mode, ...))
    }
  }
}
