function [G,U_c,V_c,W_c]    =     hosvd(Y,dim_c,type)
%-------------------------------------------------------------------------%
% [G,U_c,V_c,W_c]     =     hosvd(Y,dim_u,dim_c,type)
%
% This function approximates the HOSVD of tensor Y with approximations 
% dim_c for ranks of unfolding matrices (dimensions of the core tensor). 
% SVD can be evaluated by the standard method or by the randomized 
% approximation.
% 
% Inputs:
% 
% - Y           : data block;
% - dim_c       : dimensions of the compressed space;
% - type        : type of method to be used to evaluate the SVD.
%
%                 For type='standard', the standard method is used.
%
%                 For type='random', the randomized SVD is used.
%
% Outputs:
%
% - G           : core tensor (compressed data block);
% - U_c,V_c,W_c : HOSVD factors (compression matrices).
%
% External functions:
%
% - rsvd.m      : Randomized SVD.
%
% List of updates                 -     21/07/2015  -     J. E. Cohen
%                                       Creation of the file
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Size of the tensor
dim_u      =     size(Y);
% Index for n-mode matricization
prod_u     =     prod(dim_u)./dim_u;
%-------------------------------------------------------------------------%

%--------------------Singular Value Decompositions------------------------%
switch(type)

%---------------------------Standard SVD----------------------------------%
    case('standard')
% Standard SVD of the unfoldings to obtain the left singular vectors
% Unfolding (3)
Y_unf     =     reshape(permute(Y,[3,1,2]),dim_u(3),prod_u(3));
[W,~,~]   =     svd(Y_unf,'econ');
% Unfolding (2)
Y_unf     =     reshape(permute(Y,[2,1,3]),dim_u(2),prod_u(2));
[V,~,~]   =     svd(Y_unf,'econ');
% Unfolding (1)
Y_unf     =     reshape(Y,dim_u(1),prod_u(1));
[U,~,~]   =     svd(Y_unf,'econ');
% Truncation of the left singular vectors (low rank approximation)
U_c       =     U(:,1:dim_c(1));
V_c       =     V(:,1:dim_c(2));
W_c       =     W(:,1:dim_c(3));
%-------------------------------------------------------------------------%

%---------------------------Randomized SVD--------------------------------%
    case('random')
% RSVD of the unfoldings (one power iteration, result is already truncated)
% Unfolding (3)
Y_unf     =     reshape(permute(Y,[3,1,2]),dim_u(3),prod_u(3));
[W_c,~,~] =     rsvd(Y_unf,dim_c(3),1);
% Unfolding (2)
Y_unf     =     reshape(permute(Y,[2,1,3]),dim_u(2),prod_u(2));
[V_c,~,~] =     rsvd(Y_unf,dim_c(2),1);
% Unfolding (1)
Y_unf     =     reshape(Y,dim_u(1),prod_u(1));
[U_c,~,~] =     rsvd(Y_unf,dim_c(1),1);
%-------------------------------------------------------------------------%
otherwise
    disp('Unknown method. Choose between: standard or random.')
end

% Evaluation of the core tensor
G   =     reshape(U_c'*Y_unf,dim_c(1),dim_u(2),dim_u(3));
clear     Y_unf;
T_2 =     reshape(permute(G,[2,1,3]),dim_u(2),dim_c(1)*dim_u(3));
G   =     permute(reshape(V_c'*T_2,dim_c(2),dim_c(1),dim_u(3)),[2,1,3]);
T_3 =     reshape(permute(G,[3,1,2]),dim_u(3),dim_c(1)*dim_c(2));
G   =     permute(reshape(W_c'*T_3,dim_c(3),dim_c(1),dim_c(2)),[2,3,1]);
%-------------------------------------------------------------------------%
end

