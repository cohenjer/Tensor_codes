-------------------------M2PALS codes for computing multiple dictionary low rank factorization ----------------
------------------------------------------------------Jeremy E. Cohen and Nicolas Gillis, November 2017--------

If you use this code, please cite the following:

Spectral Unmixing with Multiple Dictionaries, Jeremy E.Cohen and Nicolas Gillis, arXiv:1711.02883

Also please do not hesitate to report any bug that you might encounter.

Installation:
Simply download the .zip folder, unzip, and add to path the following folders:
addpath .../tensor
addpath .../tools
addpath .../hierclust2nmf_v2
addpath .../nonnegative_least_squares

Some of the more involved content in the simulation folder requires the following codes:
- GLUP+SDSOMP: https://github.com/rammanouil
- MLR segmentation: http://www.umbc.edu/rssipl/people/aplaza/, under "Spectral-Spatial Hyperspectral Image Segmentation Using Subspace Multinomial Logistic Regression and Markov Random Fields"
- SNPA: https://sites.google.com/site/nicolasgillis/
- VCA: http://www.lx.it.pt/~bioucas/code.htm

Contents:
	- Main folder: 
		+ D_nmf_OMP     : computes M=AB^t where A=D(:,K)
		+ Dmult_nmf_OMP : computes M=AB^t where A=[D1(:,K1),...,Dk(:,Kk)]\Pi
		+ D_select	: a graphic interface for choosing pure pixel areas in an HSI
		+ MatchP        : finds the best K that minimizes ||B - D(:,K)||_F^2
		+ MatchP_multiD : same as MatchP but with multiple dictionaries. Uses the munkres algorithm in tools.
		+ ReadMe.txt    : this file.
	- tools folder:
		Contains functions that make developping/hacking the code easier
	- tensor folder:
		+ Dmult_mpnals  : computes the PARAFAC/CPD tensor decomposition with multiple dictionary constraints. May be buggy.
	- simulation_scripts:
		Contains codes to run the simulation of the mentioned publication. Requires some external toolboxes (see installation above).
	- nonnegative_least_squares:
		Contains functions that compute nonnegative least squares solutions in various conditions (credits to Nicolas Gillis, see his personal webpage for more information).
	- Hierclust2nmf_v2:
		Contains the code to compute the H2NMF segmentation method. Credits to Nicolas Gillis.



Copyright 2017 Jeremy Cohen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
The following reference should be cited :
Spectral Unmixing with Multiple Dictionaries, Jeremy E.Cohen and Nicolas Gillis, arXiv:1711.02883

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.