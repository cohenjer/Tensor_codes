function  [A_cor,B_cor,C_cor]     =    amb_correct(A_hat,B_hat,C_hat,A,B,C)
%-------------------------------------------------------------------------%
%         [A_cor,B_cor,C_cor]     =    amb_correct(A_hat,B_hat,C_hat,A,B,C)
% Correction of the sign and permutation ambiguity of the CP model.
%
% Inputs:
% 
% - A_hat,...   : estimated factors without correction;
% - A,...       : true factors.
%
% Outputs:
% 
% - A_cor,...   : factors with corrected sign and permutations.
%
% External functions:
%
% - munkres.m   : Hungarian algorithm.
%
% List of updates                 -     22/01/2015  -     R. C. Farias
%                                       Creation of the file
%                                 -     22/07/2015  -     R. C. Farias
%                                       Update using the Hungarian
%                                       algorithm (munkres.m)
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Dimensions
dim =     [size(A,1) size(B,1) size(C,1)];
%-------------------------------------------------------------------------%

%-------------------------Permutation correction--------------------------%
% Matrix with negative correlations between the columns
P_cost=-abs([A;B;C]'*[A_hat;B_hat;C_hat]);
% Permutation indexes given by the Hungarian algorithm
[perm,~]=munkres(P_cost);

% Permutation of the columns
A_cor     =     A_hat(:,perm);
B_cor     =     B_hat(:,perm);
C_cor     =     C_hat(:,perm);
%-------------------------------------------------------------------------%

%------------------------------Sign correction----------------------------%
%Adjusting the signs
% Factor A
sign_A    =     sign(sum(A.*A_cor));
A_cor     =     A_cor.*repmat(sign_A,dim(1),1);
 
% Factor B
sign_B    =     sign(sum(B.*B_cor));
B_cor     =     B_cor.*repmat(sign_B,dim(2),1);
 
% Factor C
sign_C    =     sign(sum(C.*C_cor));
C_cor     =     C_cor.*repmat(sign_C,dim(3),1);
%-------------------------------------------------------------------------%
end