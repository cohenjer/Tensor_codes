%---------------------------proco_als_test--------------------------------%
% This file tests the PROCO-ALS algorithm.
%
%                 For type='standard', standard nonnegativity constraints
%                 are considered.
%
%                 For type='simplex', factors are nonnegative and the
%                 factor C has a sum to one constraint on its columns.
%
%                 For type='smooth', factors are nonnegative and the factor
%                 C is smooth and can be efficiently compressed with a 
%                 bsplines basis. The bsplines basis Bs must be given.
%                 The number of nodes in the bsplines must equal dim_c(3).
%
% External functions:
%
% - kr.m              : Khatri-Rao product;
% - simplex_proj.m    : Simplex prjection algorithm;
% - proco_als.m       : PROCO-ALS algorithm;
% - amb_correct.m     : Matching Permutation for two CP models.
%
% List of updates                 -     27/07/2015  -     J. E. Cohen and
%                                                         R. C. Farias
%                                       Creation of the file
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
clc;
clear     all;
close     all;
%-------------------------------------------------------------------------%

%-----------------------------Model parameters----------------------------%
% Type of constraints
type      =     'standard';
% Uncompressed dimensions
dim_u     =     [10 10 10];
% Compressed dimensions
dim_c     =     [5 5 5];
% Number of components
R         =     5;
% Standard deviations of the measurements noise (reduce for type='simplex')
sigma_n   =     0.001;
% Number of iterations
iter      =     1000;
%-------------------------------------------------------------------------%

%------------------------------Tensor generation--------------------------%
% Unnormalized nonnegative factors
A   =     abs(randn(dim_u(1),R));
B   =     abs(randn(dim_u(2),R));
% Norm of the factors
norm_A    =     sum(A);
norm_B    =     sum(B);
% Generation of the C factor
switch(type)
    
% Standard nonnegative factor
case('standard')
C   =     abs(randn(dim_u(3),R));
% Bspline basis is set to zero
Bs  = 0;
% Renormalization of C
C         =     C.*repmat(norm_A.*norm_B,dim_u(1),1);

% Smooth nonnegative factor (bspline with n_k number of knots)
case('smooth')
% Generation of the bsplines basis
t         =     1:dim_u(3);
step      =     (t(end)-t(1))/(dim_c(3)-1);
pos       =     t(1):step:t(end);
pos_pad   =     [[pos(1)-2*step,pos(1)-step],pos,...
                [pos(end)+step,pos(end)+2*step]];
Bs        =     zeros(dim_u(3),dim_c(3));
for i=1:dim_c(3)
    Bs(:,i)     =     fnval(bspline(pos_pad(i:(i+4))),t);
    Bs(Bs<0)    =     0;
end
% Coefficients in the bsplines basis
S   =     rand(dim_c(3),R);
% C in column space of Bs
C   =     Bs*S;
% Renormalization of C
C         =     C.*repmat(norm_A.*norm_B,dim_u(1),1);

% Nonnegative C factor with sum-to-one constraint (simplex)
case('simplex')
C   =     simplex_proj(abs(rand(dim_u(3),R)));
% Bspline basis is set to zero
Bs  = 0;
otherwise
disp('Unknown constraint. Choose standard, simplex or smooth.')
end

% Normalized nonnegative factors
A         =     A./repmat(sum(A),dim_u(1),1);
B         =     B./repmat(sum(B),dim_u(2),1);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Clean and noisy data block
% Clean data block (useful to evaluate the signal to noise ratio)
X   =     reshape(A*kr(C,B)',dim_u);
% Noisy data block
Y   =     X+sigma_n*randn(dim_u);
% Signal to noise ratio
SNR =     10*log10(sum(sum(sum(X.^2)))/sum(sum(sum((Y-X).^2))));
% Delete clean data block 
clear     X;
            
%-----------------------------Initialization------------------------------%
A_0         =     abs(randn(size(A)));
B_0         =     abs(randn(size(B)));
A_0         =     A_0.*repmat(1./sum(A_0),dim_u(1),1);
B_0         =     B_0.*repmat(1./sum(B_0),dim_u(2),1);
C_0         =     simplex_proj(abs(randn(size(C))))+eps;
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%-----------------------PROCO-ALS updated version-------------------------%
tic
% ALS
[A_als,B_als,C_als,err_als_vec,U,V,W]   =     proco_als(Y,iter,dim_c,A_0,...
                                                           B_0,C_0,type,Bs);
toc

% Resulting factors
% Correcting permutation and sign ambiguities
[A_als,B_als,C_als]   =     amb_correct(A_als,B_als,C_als,A,B,C);
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

