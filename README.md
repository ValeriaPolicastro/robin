# robin <img src="man/img/logoRobin.png" align="right" height="138.5" />
 Available on CRAN https://CRAN.R-project.org/package=robin
<br/><br>

**_ROBIN (ROBustness In Network)_** is an R package for the validation of community detection. It has a double aim: it **studies the robustness** of a community detection algorithm and it **compares** the robustness of **two community detection algorithms**. 

<p align="center">
  <img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Schermata%20del%202019-09-23%2012-50-52.png" width="500" height="300" />
</p>


The package implements a methodology that detects if the community structure 
found by a detection algorithm is statistically significant or is a result 
of chance, merely due to edge positions in the network.

###### The package: 

1) **Examine the robustness** of a community detection algorithm against random perturbations of the original graph

2) **Tests the statistical difference** between the stability measure curves created

3) Makes a **comparison between different community detection algorithms** to choose the one that better fits the network of interest

4) Gives a graphical **interactive representation** 

---------------

## Example 1:   "Robustness of a community detection algorithm"
```{r}
my_network <- system.file("example/football.gml", package="robin")
graph <- prepGraph(file=my_network, file.format="gml")
graphRandom <- random(graph=graph)
proc <- robinRobust(graph=graph, graphRandom=graphRandom, measure="vi", method="louvain", type="independent")               
plotRobin(graph=graph, model1=proc$Mean, model2=proc$MeanRandom, legend=c("real data", "null model"))
```
<p align="center">
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/PlotRobin.png" width="350" height="300" />
</p>

```{r}
#For the testing:
robinFDATest(graph=graph, model1=proc$Mean, model2=proc$MeanRandom)
robinGPTest(model1=proc$Mean, model2=proc$MeanRandom)
```

## Example 2:   "Comparison of two community detection algorithms"
```{r}
my_network <- system.file("example/football.gml", package="robin")
graph <- prepGraph(file=my_network, file.format="gml")
comp <- robinCompare(graph=graph, method1="fastGreedy", method2="louvain", measure="vi", type="independent")                
plotRobin(graph=graph, model1=comp$Mean1, model2=comp$Mean2, legend=c("fastGreedy", "louvain"), title="FastGreedy vs Louvain")
```
<p align="center">
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/PlotCompare.png" width="350" height="300"/>
</p>
In this example, the Louvain algorithm fits better the network of interest, as the curve of the stability measure varies less than the one obtained by the Fast greedy method.

```{r}
#For the testing:
robinFDATest(graph=graph, model1=comp$Mean1, model2=comp$Mean2)
robinGPTest(model1=comp$Mean1, model2=comp$Mean2)
```
## Reference
ROBustness In Network (robin): an R package for Comparison and Validation of communities Valeria Policastro, Dario Righelli, Annamaria Carissimo, Luisa Cutillo, Italia De Feis. The R Journal (2021) https://journal.r-project.org/archive/2021/RJ-2021-040/index.html

## License
[Copyright (c) 2019 V. Policastro,  A. Carissimo, L. Cutillo, I. De Feis and D. Righelli.](https://github.com/ValeriaPolicastro/robin/blob/master/LICENSE)
