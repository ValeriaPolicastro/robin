# robin
ROBIN (Robustness In Network) is an R package for the Validation of community detection it has a double aim it studies the robustness of a community detection algorithm and compares the robustness of two community detection algorithms. 

The package implements a methodology that detects if the community structure 
found by a detection algorithm is statistically significant or is a result 
of chance, merely due to edge positions in the network. It performs a 
perturbation strategy and runs a null model to build a set of procedures based 
on different stability measures. 

In particular, it provides:
1)A procedure to examine the stability of the partition recovered against random 
perturbations of the original graph structure
2)Three tests to determine whether the obtained clustering departs significantly 
from the null model
3)A routine to compare different detection algorithms applied to the same 
network to discover which fits better
4)A graphical interactive representation 


Copyright (c) 2019 V. Policastro,  A. Carissimo, L. Cutillo, I. De Feis and D. Righelli.
