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
P_cost=-abs(B'*B_hat);
%P_cost=-abs([A;B;C]'*[A_hat;B_hat;C_hat]);
% Permutation indexes given by the Hungarian algorithm
[perm,~]=munkres(P_cost);

if size(B,2)<=size(B_hat,2)

% Permutation of the columns + selection
A_cor     =     A_hat(:,perm);%2 0 1 0 3
B_cor     =     B_hat(:,perm);
C_cor     =     C_hat(:,perm);

else
    
% Permutation of the columns + zero-filling
for r=1:size(B,2)
    if perm(r)==0
A_cor(:,r)     =     zeros(dim(1),1);
B_cor(:,r)     =     zeros(dim(2),1);
C_cor(:,r)     =     zeros(dim(3),1);
    else
A_cor(:,r)     =     A_hat(:,perm(r));
B_cor(:,r)     =     B_hat(:,perm(r));
C_cor(:,r)     =     C_hat(:,perm(r));
    end
end
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