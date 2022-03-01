
#' @title Visualisation function for static multiclass network models
#'
#' @param x An individual element (e.g. $linear.strs) from the output of multiJGL function.
#' @param obs.class.names Observational class names.
#' @param node.names Network node names.
#' @param obs.class.legend Optional: Obs class specific legends.
#' @param graphlayout igraph argument.
#' @param lcex the cex argument for legends.
#' @param ... Additional igraph arguments.
#'
#' @return Multiclass network plot.
#' @import jeek
#' @importFrom igraph layout_in_circle
#' @importFrom igraph subgraph.edges
#' @importFrom igraph E
#' @importFrom grDevices rainbow
#' @importFrom graphics legend
#' @importFrom igraph graph.adjacency
#' @examples print("multiJGLplot(multiJGLoutput$linear.strs)")
#' @export
multiJGLplot <- function(x, obs.class.names = NULL, node.names = NULL,
                           obs.class.legend = TRUE, graphlayout = NULL, lcex = 0.5,   ...)
{

  .env = "environment: namespace:jeek"
  obs_class.index = NULL

  adjacent.mat = returngraph(x)

  graphlayout = .makelayout(x,adjacent.mat, graphlayout = graphlayout)
  ## make title according to user input
  title = .maketitle(
    node.names = node.names, obs.class.names = NULL
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
        col = grDevices::rainbow(length(x$Graphs) + 1), pch = 18,cex = lcex)
    } else {
      legend("topleft" , legend = c(obs.class.names, "shared structures"),
             #Modify the legend position by  the lcex function argument
             col = grDevices::rainbow(length(x$Graphs) + 1), pch = 18, cex = lcex)
    }
  }
}

returngraph <- function(x, structures.to.plot = "obs.class",
                        obs_class.index = NULL, index = NULL) {
  adj = .make.adj.matrix(x$Graphs)
  diag(adj) = 0
  adjacent.mat = igraph::graph.adjacency(adj, mode = "upper", weighted = TRUE)

  K = length(x$Graphs)

  if (!is.null(igraph::E(adjacent.mat)$weight)) {
    igraph::E(adjacent.mat)$color = grDevices::rainbow(K+1)[igraph::E(adjacent.mat)$weight]
  }

  if (structures.to.plot == "share") {
    adjacent.mat = subgraph.edges(adjacent.mat, which(igraph::E(adjacent.mat)$weight == K + 1), delete.vertices = FALSE)
  } else if (structures.to.plot == "obs.class") {
    if (!is.null(obs_class.index)) {
      if (!prod(obs_class.index %in% (0:K))) {
        stop("please specify valid obs.class number(s)")
      }
      adjacent.mat = subgraph.edges(adjacent.mat, which(igraph::E(adjacent.mat)$weight %in% c(obs_class.index, K + 1)), delete.vertices = FALSE)
    }
  }

  return(adjacent.mat)
}


### helper function to make title (see also jeek R package).
.maketitle <- function(structures.to.plot = "obs.class", obs_class.index = NULL,
                       index = NULL, node.names = NULL, obs.class.names = NULL)
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
      return (paste("Estimated", "network structures", sep = " "))
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


