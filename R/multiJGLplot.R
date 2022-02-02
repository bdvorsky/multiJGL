
#' @title Visualisation function for static multiclass network models
#'
#' @param x Estimated network structures from the static_LNLJGL function.
#' @param obs.class.names Observational class names.
#' @param node.names Network node names.
#' @param network.type Should linear, nonlinear or combined networks be plotted with options: "linear", "nonlinear" or "both"
#' @param structures.to.plot Should obs.class specific structures to be specified.
#' @param obs.class.legend Optional: Obs class specific legends.
#' @param obs_class.index  Optional: Used to specify the subnetworks.
#' @param graphlayout igraph argument.
#' @param lcex the cex argument for legends.
#' @param ... Additional igraph arguments.
#'
#' @return Multiclass network plot.
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @import jeek
#' @import igraph
#' @export
#'
#' @examples print("multiJGLplot(x)")
NULL
utils::globalVariables(c("network.type"))



multiJGLplot <- function(x, obs.class.names = NULL, node.names = NULL,
                         network.type = "both",
                         structures.to.plot = "obs.class",
                         obs.class.legend = TRUE,
                         obs_class.index = NULL, graphlayout = NULL, lcex = 0.5,   ...)
{

  .env = "environment: namespace:jeek"
  obs_class.index = unique(obs_class.index)


  adjacent.mat = returngraph(x, structures.to.plot = structures.to.plot,

                             obs_class.index = obs_class.index
  )

    graphlayout = .makelayout(x,adjacent.mat, graphlayout = graphlayout)
  ## make title according to user input

  title = .maketitle(
    structures.to.plot = structures.to.plot,
    obs_class.index = obs_class.index,
    node.names = node.names,  nettype = network.type, obs.class.names = NULL
  )

  plot(
    adjacent.mat,
    layout = graphlayout,
    vertex.label.font = 2,
    vertex.shape = "none",
    vertex.label.color = "gray40",
    vertex.label = node.names,
    vertex.label.cex = .7,
    vertex.frame.color = "white",
    vertex.size = 10 ,
    main = title,
    ...
  )

  if(obs.class.legend == TRUE){

    if(is.null(obs.class.names)){
      legend("topleft" , legend = c(paste("obs.class", c(
        1:length(x$Graphs)), "specific"), "shared structures"),
        col = rainbow(length(x$Graphs) + 1), pch = 18,cex = lcex)
    } else {
      legend("topleft" , legend = c(obs.class.names, "shared structures"),
             #Modify the legend position by  the lcex function argument
             col = rainbow(length(x$Graphs) + 1), pch = 18, cex = lcex)
    }
  }
}

returngraph <- function(x, structures.to.plot = "obs.class",
                        obs_class.index = NULL, index = NULL) {
  adj = .make.adj.matrix(x$Graphs)
  diag(adj) = 0
  adjacent.mat = graph.adjacency(adj, mode = "upper", weighted = TRUE)

  K = length(x$Graphs)

  if (!is.null(E(adjacent.mat)$weight)) {
    E(adjacent.mat)$color = rainbow(K+1)[E(adjacent.mat)$weight]
  }

  if (structures.to.plot == "share") {
    adjacent.mat = subgraph.edges(adjacent.mat, which(E(adjacent.mat)$weight == K + 1), delete.vertices = FALSE)
  } else if (structures.to.plot == "obs.class") {
    if (!is.null(obs_class.index)) {
      if (!prod(obs_class.index %in% (0:K))) {
        stop("please specify valid obs.class number(s)")
      }
      adjacent.mat = subgraph.edges(adjacent.mat, which(E(adjacent.mat)$weight %in% c(obs_class.index, K + 1)), delete.vertices = FALSE)
    }
  }

  return(adjacent.mat)
}


### helper function to make title (see also jeek R package).
.maketitle <- function(structures.to.plot = "obs.class", obs_class.index = NULL,
                       index = NULL,  nettype = network.type, node.names = NULL, obs.class.names = NULL)
{
  if (structures.to.plot == "share") {
    return ("shared structures")
  }

  if (structures.to.plot == "obs.class_specific") {
    temp = paste(as.character(obs_class.index), collapse = ", ")
    return (paste("Obs.class", temp, "specific network"))
  }

  if (structures.to.plot == "obs.class") {
    if (is.null(obs_class.index)) {
      return (paste("Estimated", nettype, "network structures", sep = " "))
    }
    else {
      if (length(obs_class.index) == 1) {
        if (obs_class.index == 0) {
          return ("Common network structures")
        } else {
          return (paste("obs.class", obs_class.index, "network"))
        }
      } else {
        if (0 %in% obs_class.index) {
          temp = obs_class.index[-(which(obs_class.index %in% 0))]

          return(paste("obs.class", paste(as.character(temp), collapse = ", ")), "network")
        } else {
          return (paste("obs.class", paste( as.character(obs_class.index), collapse = ","), "network"))
        }
      }
    }
  }
}

### helper function to create layout for graph
.makelayout <-
  function(x,adjacent.mat,  graphlayout = NULL){
    if (is.null(graphlayout)) {
      graphlayout = layout_in_circle(adjacent.mat, order = rev(c(1:dim(x$Graphs[[1]])[1])))
    }
    return(graphlayout)
  }

### a jeek R package helper function to make adjacent matrix
.make.adj.matrix <-  function(theta, separate=FALSE){
  K = length(theta)
  adj = list()
  if(separate) {
    for(k in 1:K){
      adj[[k]] = (abs(theta[[k]])>1e-5)*1
    }
  }
  if(!separate){
    adj = 0*theta[[1]]

    for(k in 1:K){
      adj = adj+(abs(theta[[k]])>1e-5)*2^(k-1)
    }
  }
  return(adj)
}





