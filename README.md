# robin
ROBIN (Robustness In Network) gives a statistical answer to the validation of 
the community structure by looking at the robustness of the network and 
compares the robustness of two community detection algorithms to see which fits 
better the network of interest.

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
