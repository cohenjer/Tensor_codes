# Tensor_codes
MATLAB codes for computing various tensor decomposition. Most of the shared code is rather unoptimized, to be used to check viability of the proposed new tensor decompositions models. Most algorithms are also based on variants of Alternating Least Squares.

Last Update : 05/04/2018

-------- Contents ---------
---------------------------

1/ Dictionary-based decompositions
----------------------------------
  A set of functions to decompose a tensor using a CPD model, where one factor lives among a large dictionary of known components.
  
  a) M2PALS : Multiple dictionary are available, with bounds on the number of atoms to select for each dictionary.
  
  b) MPALS : factor A in the CPD of a tensor T is so that A = D(:,K), K a set of idexes. Features greedy and flexible algorithms.
  
2/ Coupled Decompositions
-------------------------

  a) CCP : Flexibly coupled tensor decompositions. 
  
  b) NNP2 : Flexibly coupled PARAFAC2 with nonnegativity constraints in the coupled mode. 
  
3/ Constrained Compression and acceleration for constrained tensor data
-----------------------------------------------------------------------

  a) PROCO-ALS : Fast nonnegative tensor PARAFAC/Canonical Polyadic decomposition. Compression is based on randomized SVD.
  
4/ Nonlinear tensor decomposition 
---------------------------------

  a) NLFD : Nonlinear fluorescence decomposition, designed for fluorescence samples with high concentrations, where the linearity of the CPD model does not hold. Based on the Levenberg Marquardt algorithm.
