----- READ ME-----

Main files :
- comparison_test_file_ICASSP_final : reproduces Figure 2 of the paper (runs for a few hours)
- comparison_test_file_ICASSP_revision : reproduces the figures of the response to reviewers (fast)

To run the experiments of the ICASSP paper, you need the following additional toolboxes
- ksubspace : https://fr.mathworks.com/matlabcentral/fileexchange/37353-k-subspaces?focused=5236203&tab=function
- sNNOMP: ask the authors of `` An optimized version of non-negative OMP'', GRETSI 2017

To allow the experiments to run without these codes, some lines of code are commented in the scripts. 
Please uncomment the required lines if these toolboxes are in your path.

Also, the following codes are provided as side functions:
- nnlsHALSupdt, nnlsHALSupdtv2, ... --> NNLS solvers, with various constraints (l1, fixed support, l2...)
- HALSacc --> computes NMF using nnlsHALSupdt
- sHALSreview ---> NMF with l1 constraints
- as_k_s_nmf + solveNormalEqComb ---> active-set sparse NMF (external code)
- k_summits ---> proposed NOLRAK algorithm
- k_s_nmf ---> exact combinatorial sparse NMF (ESNA)

The code is provided as is, without any guaranty. If you use either k_s_nmf or k_summits provided here, please cite:
"Jeremy E. Cohen, Nicolas Gillis, Nonnegative Low-rank Sparse Component Analysis, ICASSP2019"

