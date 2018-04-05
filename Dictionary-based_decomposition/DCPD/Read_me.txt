-------------------------MPALS codes for computing dictionary-based low rank factorization --------------------
------------------------------------------------------Jeremy E. Cohen and Nicolas Gillis, November 2017--------

If you use this code, please cite the following:

"Dictionary-based Tensor Canonical Polyadic Decomposition, J.E.Cohen, N.Gillis, IEEE Transactions on Signal Processing, vol.66 issue 7, pp. 1876-1889 2018"

Also please do not hesitate to report any bug that you might encounter.

Installation:
Simply download the DCPD.zip folder, unzip, and add to path the following folders:
addpath .../Continuous
addpath .../Simulations
addpath .../Greedy
addpath .../Tools
addpath .../Matrix

Some of the more involved content in the simulation folder requires the following codes:
- GLUP+SDSOMP: https://github.com/rammanouil
- SNPA: https://sites.google.com/site/nicolasgillis/
- NFINDR: For instance, http://davidkun.github.io/

Contents:
	- Main folder: 
		+ D_mpnals        : computes T = I x A x B x C  where B=D(:,K) using ALS and a greedy projection technique (can be used without nonnegative constraints)
		+ D_nmf_OMP 	  : similar function but for matrix data M=AB^t (two way arrays) where A=D(:,K)
		+ D_mpnals_flex_2 : computes T = I x A x B x C  where B \approx D(:,K), parameters inc is delta := delta*inc as long as ||B-D(:,K)||_F/\|B\|_F > lambda. Should converge when delta is fixed.
		+ D_nmf_FlexMP    : similar function but for matrix data. Here, inc is set to 1.5 by default and lambda is called beta.

	- TOOLS folder:
		Contains functions that make developping/hacking the code easier. 
	- Greedy folder:
		Contains Smooth-MPALS for high order tensor data, and not fully developped functions for Dictionary Parafac2 CPD
	- Continuous folder:
		Contains the fast-gradient algorithm for high order tensor data. Not recommanded for practical use.
	- Matrix folder:
		Contains Smooth-MPALS for matrix data, as well as run files for computing specral unmixing of Urban, Terrain and San Diego Airport Data set.
		**WARNING: To run these simulations, unzip the SU.Data_CodesGillis file and add to path folders and subfolders. It contains the data and some comparison methods like SPA and SNPA (credits to Nicolas Gillis).**
	- Simulations:
		Contains the code to reproduce the simulations of the article "Dictionary-based Tensor Canonical Polyadic Decomposition".

Copyright 2017 Jeremy Cohen

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
The following reference should be cited :
"Dictionary-based Tensor Canonical Polyadic Decomposition, J.E.Cohen, N.Gillis, IEEE Transactions on Signal Processing, vol.66 issue 7, pp. 1876-1889 2018"

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
