function  [A_cor,B_cor,C_cor]     =    amb_correct_1f(A_p,B_p,C_p,C)
%-------------------------------------------------------------------------%
%         [A_cor,B_cor,C_cor]  =    amb_correct_1f(A_hat,B_hat,C_hat,A,B,C)
% Correction of the sign and permutation ambiguity for the third factor
% of estimated coupled CP models.
%
% Inputs:
% 
% - A_p,B_p,C_p : estimated factors without correction;
% - C           : estimated coupled factor of the first tensor.
%
% Outputs:
% 
% - A_cor,...    : factors with corrected sign and permutations.
%
% External functions:
%
% - munkres.m   : Hungarian algorithm.
%
% List of updates                 -     23/02/2015  -     Rodrigo Cabral
%                                       Creation of the file
%                                 -     22/07/2015  -     R. C. Farias
%                                       Update using the Hungarian
%                                       algorithm (munkres.m)
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Dimensions
dim =     [size(A_p,1) size(B_p,1) size(C_p,1)];
%-------------------------------------------------------------------------%

%-------------------------Permutation correction--------------------------%
% Matrix with negative correlations between the columns
P_cost=-abs(C'*C_p);
% Permutation indexes given by the Hungarian algorithm
[perm,~]=munkres(P_cost);

% Permutation of the columns
A_cor     =     A_p(:,perm);
B_cor_init=     B_p(:,perm);
C_cor_init=     C_p(:,perm);
%-------------------------------------------------------------------------%

%------------------------------Sign correction----------------------------%
%Adjusting the signs
% Factor C
sign_C    =     sign(sum(C.*C_cor_init));
C_cor     =     C_cor_init.*repmat(sign_C,dim(3),1);
B_cor     =     B_cor_init.*repmat(sign_C,dim(2),1);
%-------------------------------------------------------------------------%
end