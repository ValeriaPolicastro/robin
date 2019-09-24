# robin <img align= "right" src="https://github.com/ValeriaPolicastro/robin/blob/master/images/Schermata%20del%202019-09-23%2016-15-54.png" width="100" height="100" /> 

**_ROBIN (Robustness In Network)_** is an R package for the Validation of community detection it has a double aim it studies the robustness of a community detection algorithm and compares the robustness of two community detection algorithms. 

<p align="center">
  <img src="https://github.com/ValeriaPolicastro/robin/blob/master/images/Schermata%20del%202019-09-23%2012-50-52.png" width="500" height="300" />
</p>


The package implements a methodology that detects if the community structure 
found by a detection algorithm is statistically significant or is a result 
of chance, merely due to edge positions in the network.

###### It provides:
* A procedure to examine the stability of the partition recovered against random 
perturbations of the original graph structure
* Three tests to determine whether the obtained clustering departs significantly 
from the null model
* A routine to compare different detection algorithms applied to the same 
network to discover which fits better
* A graphical interactive representation
---------------

## Example 1: "Robustness of a community detection"
```{r}
my_network <- system.file("example/football.gml", package="robin")
graph <- prepGraph(file=my_network, file.format="gml")
graphRandom <- random(graph=graph)
proc <- robinRobust(graph=graph, graphRandom=graphRandom, measure="vi", 
                  method="louvain", type="independent")               
plotRobin(graph=graph, model1=proc$Mean, model2=proc$MeanRandom, 
legend=c("real data", "null model"), measure="vi")
robinGPTest(ratio=proc$ratios)
```
<p align="center">
<img src="https://github.com/ValeriaPolicastro/robin/blob/master/Figures%20Paper/PlotRobin.jpeg" width="400" height="250" />
</p>

## Example 2: "Comparison of two community detection"
```{r}
my_network <- system.file("example/football.gml", package="robin")
graph <- prepGraph(file=my_network, file.format="gml")
comp <- robinCompare(graph=graph, method1="fastGreedy",
                method2="louvain", measure="vi", type="independent")                
plotRobin(graph=graph, model1=comp$Mean1, model2=comp$Mean2, measure="vi", 
legend=c("fastGreedy", "louvain"), title="FastGreedy vs Louvain")
robinAUC(graph=graph, model1=comp$Mean1, model2=comp$Mean2, measure="vi")
```
<p align="center">
<img src="https://github.com/ValeriaPolicastro/robin/blob/master/Figures%20Paper/PlotCompare.jpeg" width="400" height="250"/>
</p>

## License
[Copyright (c) 2019 V. Policastro,  A. Carissimo, L. Cutillo, I. De Feis and D. Righelli.](https://github.com/ValeriaPolicastro/robin/blob/master/LICENCE)

