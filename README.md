<<<<<<< HEAD
# robin <img align= "right" src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Logo2.png" width="70" height="70" /> 

**_ROBIN (ROBustness In Network)_** is an R package for the validation of community detection it has a double aim it **studies the robustness** of a community detection algorithm and **compares** the robustness of **two community detection algorithms**. 

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
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Schermata%20del%202019-09-23%2012-24-29.png" width="500" height="350" />
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
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Schermata%20del%202019-09-23%2012-34-23.png" width="500" height="350"/>
</p>
In this example, the Louvain algorithm fits better the network of interest, as the curve of the stability measure varies less than the one obtained by the Fast greedy method.

## License
[Copyright (c) 2019 V. Policastro,  A. Carissimo, L. Cutillo, I. De Feis and D. Righelli.](https://github.com/ValeriaPolicastro/robin/blob/master/LICENSE)

=======
# robin <img align= "right" src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Logo2.png" width="70" height="70" /> 

**_ROBIN (ROBustness In Network)_** is an R package for the validation of community detection it has a double aim it **studies the robustness** of a community detection algorithm and **compares** the robustness of **two community detection algorithms**. 

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
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Schermata%20del%202019-09-23%2012-24-29.png" width="500" height="350" />
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
<img src="https://github.com/ValeriaPolicastro/Paper-Robin/blob/master/images/Schermata%20del%202019-09-23%2012-34-23.png" width="500" height="350"/>
</p>
In this example, the Louvain algorithm fits better the network of interest, as the curve of the stability measure varies less than the one obtained by the Fast greedy method.

## License
[Copyright (c) 2019 V. Policastro,  A. Carissimo, L. Cutillo, I. De Feis and D. Righelli.](https://github.com/ValeriaPolicastro/robin/blob/master/LICENSE)

>>>>>>> 7846e4233edb45c61ddf2d43e1ac9d1255d22411
