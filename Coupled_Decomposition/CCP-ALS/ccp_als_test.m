%------------------------------ccp_als_test-------------------------------%
% This file tests the Coupled CP ALS algorithms.
%
% External functions:
%
% - kr.m              : Khatri-Rao product;
% - als.m             : ALS algorithm;
% - amb_correct_1f.m  : Matching Permutation for two CP models based on 
%                       factors C only;
% - amb_correct.m     : Matching Permutation for two CP models;
% - ccp_als.m         : Coupled CP ALS.
%
% List of updates                 -     22/07/2015  -     R. C. Farias 
%                                       Creation of the file
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
clc;
clear     all;
close     all;
%-------------------------------------------------------------------------%

%-----------------------------Model parameters----------------------------%
% 3-way block dimensions
dim       =     [50 50 50];
% Number of components
R         =     3;
% Standard deviations of the measurements noise
sigma_n   =     [0.2 0.001];
% Standard deviation of the coupling
sigma_c   =     0.001;
% Number of iterations
iter      =     100;
%-------------------------------------------------------------------------%

%------------------------------Tensor generation--------------------------%
% Random factors generation
A         =     randn(dim(1),R);
B         =     randn(dim(2),R);
A_p       =     randn(dim(1),R);
B_p       =     randn(dim(2),R);
C_p       =     randn(dim(3),R);
% Shared factor
C         =     C_p;
% Normalization
A         =     A.*repmat(1./sqrt(sum(A.^2)),dim(1),1);
B         =     B.*repmat(1./sqrt(sum(B.^2)),dim(2),1);
A_p       =     A_p.*repmat(1./sqrt(sum(A_p.^2)),dim(1),1);
B_p       =     B_p.*repmat(1./sqrt(sum(B_p.^2)),dim(2),1);
% Noise on the coupling
% Shared factor
C         =     C+sigma_c*randn(dim(3),R);
%-------------------------------------------------------------------------%

%---------------------------------Coupling model--------------------------%
% Coupling matrices
H   =     eye(dim(3));
H_p =     eye(dim(3));
% Coupling matrices for the vectorized form
H_vec     =     kron(eye(R),H);
H_p_vec   =     kron(eye(R),H_p);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
%Noisy measurements
% First noisy tensor
Y_1       =     reshape(A*transpose(kr(C,B)),dim)+...
                sigma_n(1)*randn(dim);
% Second noisy tensor
Y_2       =     reshape(A_p*transpose(kr(C_p,B_p)),dim)+...
                sigma_n(2)*randn(dim);

%-----------------------------Initialization------------------------------%
% First tensor
A_0         =     randn(size(A));
B_0         =     randn(size(B));
A_0         =     A_0.*repmat(1./sqrt(sum(A_0.^2)),dim(1),1);
B_0         =     B_0.*repmat(1./sqrt(sum(B_0.^2)),dim(2),1);
% Second tensor   
A_p_0       =     randn(size(A_p));
B_p_0       =     randn(size(B_p));
A_p_0       =     A_p_0.*repmat(1./sqrt(sum(A_p_0.^2)),dim(1),1);
B_p_0       =     B_p_0.*repmat(1./sqrt(sum(B_p_0.^2)),dim(2),1);
 
C_p_0       =     randn(size(C_p));
C_0         =     randn(size(C));
%-------------------------------------------------------------------------%

%---------------------------Uncoupled ALS---------------------------------%
% Uncoupled ALS for the first tensor
[A_u,B_u,C_u,err_vec_1]     =     als(Y_1,iter,A_0,B_0,C_0);
% Uncoupled ALS for the second tensor
[A_p_u,B_p_u,C_p_u,err_vec_2]     =     ...
                                        als(Y_2,iter,A_p_0,B_p_0,C_p_0);
                                    
% Correcting permutation for the initialization of the coupled algorithms
[A_p_u,B_p_u,C_p_u]   =           amb_correct_1f(A_p_u,B_p_u,C_p_u,C_u);

% Resulting factors
% Correcting permutation and sign ambiguity for the first tensor
[A_u,B_u,C_u]         =     amb_correct(A_u,B_u,C_u,A,B,C);
% Correcting permutation and sign ambiguity for the second tensor
[A_p_u,B_p_u,C_p_u]   =     ...
                            amb_correct(A_p_u,B_p_u,C_p_u,A_p,B_p,C_p);
%-------------------------------------------------------------------------%

%--------------Initialization with uncoupled results----------------------%
% First CP model
A_0 =     A_u;
B_0 =     B_u;
C_0 =     C_u;
% Second CP model
A_p_0     =     A_p_u;
B_p_0     =     B_p_u;
C_p_0     =     C_p_u;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%----------------CCP ALS updated version (Joint version)------------------%
tic
% CCP_ALS
[A_c_j,B_c_j,C_c_j,A_p_c_j,B_p_c_j,C_p_c_j,map_j_vec]  = ...
 ccp_als(Y_1,Y_2,iter,H_vec,H_p_vec,sigma_n,sigma_c,'joint',A_0,B_0,C_0,...
                                                        A_p_0,B_p_0,C_p_0);
toc
% Resulting factors
% Correcting permutation and sign ambiguity for the first tensor
[A_c_j,B_c_j,C_c_j]   =     ...
                            amb_correct(A_c_j,B_c_j,C_c_j,A,B,C);
% Correcting permutation and sign ambiguity for the second tensor
[A_p_c_j,B_p_c_j,C_p_c_j]   =     amb_correct(A_p_c_j,...
                                        B_p_c_j,C_p_c_j,A_p,B_p,C_p);
%-------------------------------------------------------------------------%

% %------------CCP ALS updated version (Alternating version)----------------%
% tic
% % CCP_ALS
% [A_c_a,B_c_a,C_c_a,A_p_c_a,B_p_c_a,C_p_c_a,map_a_vec]  = ...
%  ccp_als(Y_1,Y_2,iter,H,H_p,sigma_n,sigma_c,'alter',A_0,B_0,C_0,...
%                                                         A_p_0,B_p_0,C_p_0);
% toc
% % Resulting factors
% % Correcting permutation and sign ambiguity for the first tensor
% [A_c_a,B_c_a,C_c_a]   =     ...
%                             amb_correct(A_c_a,B_c_a,C_c_a,A,B,C);
% % Correcting permutation and sign ambiguity for the second tensor
% [A_p_c_a,B_p_c_a,C_p_c_a]   =     amb_correct(A_p_c_a,...
%                                         B_p_c_a,C_p_c_a,A_p,B_p,C_p);
% %-------------------------------------------------------------------------%

%------------CCP ALS updated version (Direct coupling version)------------%
tic
% CCP_ALS
[A_c_d,B_c_d,C_c_d,A_p_c_d,B_p_c_d,C_p_c_d,map_d_vec]  = ...
 ccp_als(Y_1,Y_2,iter,H,H_p,sigma_n,sigma_c,'direct',A_0,B_0,C_0,...
                                                        A_p_0,B_p_0,C_p_0);
toc
% Resulting factors
% Correcting permutation and sign ambiguity for the first tensor
[A_c_d,B_c_d,C_c_d]   =     ...
                            amb_correct(A_c_d,B_c_d,C_c_d,A,B,C);
% Correcting permutation and sign ambiguity for the second tensor
[A_p_c_d,B_p_c_d,C_p_c_d]   =     amb_correct(A_p_c_d,...
                                        B_p_c_d,C_p_c_d,A_p,B_p,C_p);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

