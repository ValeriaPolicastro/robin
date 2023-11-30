
################ PLOT GRAPH ###############
#' plotGraph
#'
#' @description Graphical interactive representation of the network.
#' @param graph The output of prepGraph.
#'
#' @return Creates an interactive plot, a D3 JavaScript network graph.
#' @import networkD3
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' plotGraph (graph)
plotGraph <- function(graph)
{
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph)
    plot <- networkD3::simpleNetwork(graph_d3$links, opacity=0.8, zoom=TRUE,
                                     nodeColour = "#2E66AC",
                                     fontSize=12)
    return(plot)
}   


######################## PLOT COMMUNITIES ##############
#' plotComm
#' 
#' @description Graphical interactive representation of the network and its 
#' communities.
#' 
#' @param graph The output of prepGraph.
#' @param members A membership vector of the community structure, the output of
#' membershipCommunities. 
#'
#' @return Creates an interactive plot with colorful communities, a D3 
#' JavaScript network graph.
#' @import networkD3 
#' @importFrom methods is
#' @export
#'
#' @examples
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' members <- membershipCommunities (graph=graph, method="louvain")
#' plotComm(graph, members)
plotComm <- function(graph, members) 
{
    stopifnot(methods::is(graph, "igraph"))
    stopifnot(methods::is(members, "membership"))
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph, group=members)
    # Create network plot
    plot <- networkD3::forceNetwork(Links=graph_d3$links, Nodes=graph_d3$nodes,
                                    Source='source', Target='target', 
                                    NodeID='name', Group='group', opacity=0.8,
                                    fontSize=12, legend=TRUE)
    return(plot)
}


############PLOT##############
#' plot.robin
#'
#' @description This function plots two curves: the measure of the null model and the measure
#' of the real graph or the measure of the two community detection algorithms.
#' @param x A robin class object. The output of the functions:  
#' \code{\link{robinRobust}} and \code{\link{robinCompare}}.
#' @param title The title for the graph. The default is "Robin plot".
#' @param ... other parameter
#'
#' @return A ggplot object.
#' @import ggplot2 igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' comp <- robinCompare(graph=graph, method1="fastGreedy",method2="infomap")
#' plot(comp)
#' 
plot.robin <- function(x, title="Robin plot", ...)
{   
    
    legend <- c(x$model1,x$model2)
    if (length(x$Mean1)==0)
    {
        model1 <- x$Mean
        model2 <- x$MeanRandom 
    }else{
         model1 <- x$Mean1
         model2 <- x$Mean2 
    }
   
    mvimodel1 <- as.vector((apply(model1, 2, mean)))
    mvimodel2 <- as.vector((apply(model2, 2, mean)))
    
    
    percPert <- rep((seq(0,60,5)/100), 2)
    mvi <- c(mvimodel1, mvimodel2)
    model <-c(rep("model1",13),rep("model2",13))
    dataFrame <- data.frame(percPert, mvi, model)
    plotModel <- ggplot2::ggplot(dataFrame, aes(x = percPert, 
                                                y = as.numeric(as.character(mvi)),
                                                colour = model, 
                                                group = factor(model))) +
        geom_line() +
        geom_point() +
        xlab("Percentage of perturbation") +
        ylab("Measure") +
        ggplot2::ylim(0,1)+
        ggtitle(title)
    
    cols <- c("model1" = "#00BFC4", "model2" = "#F8766D")
    plot <-  plotModel+ggplot2::scale_colour_manual(values = cols,
                                                    breaks = c("model1", "model2"), 
                                                    labels=c(legend[1], legend[2]))
    
    return(plot)
}