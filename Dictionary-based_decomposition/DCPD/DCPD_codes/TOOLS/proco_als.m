function [A,B,C,err_als_vec,U,V,W]    =     proco_als(Y,iter,dim_c,type,...
                                                            A_0,B_0,C_0,Bs)
%-------------------------------------------------------------------------%
% [A_c,B_c,C_c,err_als_vec]=     proco_als(Y,iter,dim_c,type,A_0,B_0,C_0)
%
% Compressed and projected ALS algorithm for large-scale 3-way nonnegative
% tensor decomposition.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - dim_c       : vector with the dimensions of the compressed model;
% - type        : type of constraint in the uncompressed domain.
%
%                 For type='standard', standard nonnegativity constraints
%                 are considered.
%
%                 For type='simplex', factors are nonnegative and the
%                 factor C has a sum to one constraint on its columns.
%
%                 For type='smooth', factors are nonnegative and the factor
%                 C is smooth and can be efficiently compressed with a 
%                 B-splines basis. The B-splines basis Bs must be given.
%                 The number of nodes in the B-splines must equal dim_c(3).
% 
% - A_0,B_0,C_0 : initial parameters for the factors.
%
% Outputs:
% 
% - A           : first factor;
% - B           : second factor;
% - C           : third factor;
% - err_als_vec : vector of reconstruction error on the compressed tensor.
%
% External functions:
%
% - kr.m              : Khatri-Rao product.
% - hosvd.m           : Higher-order SVD.
% - simplex_proj.m    : Projection on the simplex.
%
% List of updates                 -     28/08/2014  -     R. C. Farias
%                                       Creation of the file
%                                 -     02/10/2014  -     J. E. Cohen
%                                       Reshaping, modified outputs
%                                 -     06/07/2015  -     J. E. Cohen
%                                       Some modifications
%                                 -     23/07/2015  -     R. C. Farias
%                                       Added hosvd function and faster
%                                       computation of the error
%                                 -     11/08/2015  -     J. E. Cohen
%                                       Normalisation in the uncompressed
%                                       space and modification of simplex
%                                       constraint normalisation.
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     1;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     eps;%sqrt(eps);
stop_err  =     eps;%sqrt(eps);
% Size of the tensor
dim_u     =     size(Y);
% Index for n-mode matricization
prod_c    =     prod(dim_c)./dim_c;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
%-------------------------------------------------------------------------%

%-------------------Compression and initialization------------------------%
% Truncated HOSVD
[G,U,V,W] =     hosvd(Y,dim_c,'random');

if  strcmp(type,'smooth')
% Index for n-mode matricization
prod_u     =     prod(dim_u)./dim_u;
% Compression factor from B-splines basis
[W,Rq]    =     qr(Bs);
W         =     W(:,1:dim_c(3));
% Matrix of coefficients on the B-splines basis
Rq        =     Rq(1:dim_c(3),:);
invRq     =     pinv(Rq);

% Evaluation of the core tensor
G   =     reshape(U'*reshape(Y,dim_u(1),prod_u(1)),dim_c(1),dim_u(2),...
          dim_u(3));
T_2 =     reshape(permute(G,[2,1,3]),dim_u(2),dim_c(1)*dim_u(3));
G   =     permute(reshape(V'*T_2,dim_c(2),dim_c(1),dim_u(3)),[2,1,3]);
T_3 =     reshape(permute(G,[3,1,2]),dim_u(3),dim_c(1)*dim_c(2));
G   =     permute(reshape(W'*T_3,dim_c(3),dim_c(1),dim_c(2)),[2,3,1]);
end

% Three different unfoldings of the compressed tensor
G_mode    =     cell(3);
% Squared norm of G
norm_G          =     sqrt(G(:)'*G(:));
% Normalized G
G               =     G./norm_G;
% Data tensor unfolding
G_mode{1} =     reshape(permute(G,[1,2,3]),dim_c(1),prod_c(1));  
G_mode{2} =     reshape(permute(G,[2,1,3]),dim_c(2),prod_c(2));
G_mode{3} =     reshape(permute(G,[3,1,2]),dim_c(3),prod_c(3));
% Storing initial condition on the compressed space
A_c =    (U')*A_0; 
B_c =    (V')*B_0;
C_c =    (W')*C_0;
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     G_mode{1}-(A_c)*transpose(kr(C_c,B_c));
err_als_old     =     (res(:)')*res(:);
% Vector of squared norms of the error
err_als_vec     =     zeros(1,iter+1);
% First error
err_als_vec(1)  =     err_als_old;
fprintf('\nPROCO-ALS\n')
fprintf('----\n')
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%

%--------------------------Algorithm iterate------------------------------%
t         =     0;
while     (mean(diff_err_filt)>stop_diff &&  t<iter ...
          &&  err_als_old>stop_err   )
    % Iteration counter (counts/iter_err)
    t=    t+1;
    % PROCO-ALS (iter_err iterates)
    for   i           =     1:iter_err
    
    % PROCO solution for each factor

    % A factor update
    temp_CB     =  (C_c'*C_c).*(B_c'*B_c);
    A_c         =  (G_mode{1}*kr(C_c,B_c))/temp_CB;
    % Projection to have nonnegative factors on the uncompressed space
    A_u   =     U*A_c;          % Uncompressed estimated factor 
    A_u(A_u<0)  =     eps;        % Projection on the nonnegative orthant     
    A_u         =     A_u.*repmat(1./sum(A_u),dim_u(1),1);
    A_c         =     U'*A_u;   % Return compressed estimated factor

    
    % B factor update
    temp_CA     =  (C_c'*C_c).*(A_c'*A_c);
    B_c         =  (G_mode{2}*kr(C_c,A_c))/temp_CA;
    % Projection to have nonnegative factors on the uncompressed space
    B_u   =     V*B_c;
    B_u(B_u<0)  =     eps;
    %sum(B_u)
    if  strcmp(type,'simplex')==0
    B_u         =     B_u.*repmat(1./sum(B_u),dim_u(2),1);
    end
    B_c         =     V'*B_u;
        
    % C factor update
    temp_BA     =     (B_c'*B_c).*(A_c'*A_c);
    G_BA        =     G_mode{3}*kr(B_c,A_c);
    C_c         =     G_BA/temp_BA;
    % Projection on the uncompressed space
    switch(type)
        
    % Standard PROCO-ALS (nonnegative) 
    case('standard')
    C_u   =     W*C_c;
    C_u(C_u<0)  =     0;
    C_c         =     W'*C_u;
    
    % PROCO-ALS with sum to one constraint on C 
    case('simplex')
    C_u   =     W*C_c;
    C_u   =     simplex_proj(C_u);
    C_c   =     W'*C_u;
    
    % PROCO-ALS with smooth columns on C 
    case('smooth')
    C_u   =     invRq*C_c;
    C_u(C_u<0)  =     0;
    C_c         =     Rq*C_u;
    otherwise
    disp('Unknown constraint. Choose standard, simplex or smooth.')
    end
    end

    % Squared norm of the error (it uses some precalculated quantities)
    err_als           =     1-2*sum(sum(G_BA.*C_c))+...
                            sum(sum((C_c'*C_c).*temp_BA));
    
    % Evaluation of the relative error decrease
    diff_err_filt(1:end-1)  =     diff_err_filt(2:end);
    diff_err_filt(end)      =     abs((err_als-err_als_old)/err_als);
    
    % Storing the last error
    err_als_old       =     err_als;
    err_als_vec(t+1)  =     err_als;
    
    if mod(t*iter_err,100)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_als,t*iter_err)
    end
end

% Uncompressed factors A and B
A         =     U*A_c;
B         =     V*B_c;
C         =     W*C_c;
% Projection
A(A<0)    =     0;
B(B<0)    =     0;  
% Normalization
norm_A_u  =     sum(A);
norm_B_u  =     sum(B);
% Uncompressed and projected C factor
if strcmp(type,'simplex')
A   =     A.*repmat(norm_G*norm_B_u,dim_u(1),1);
B   =     B.*repmat(1./norm_B_u,dim_u(2),1);
C   =     simplex_proj(C);
else

A   =     A.*repmat(1./norm_A_u,dim_u(1),1);
B   =     B.*repmat(1./norm_B_u,dim_u(2),1);

C         =     W*C_c;
C(C<0)    =     0;
% Normalization
C   =     C.*repmat(norm_G*norm_A_u.*norm_B_u,dim_u(3),1);
end
% Dealing with NaN
A(isnan(A))     =     0;
B(isnan(B))     =     0;
C(isnan(C))     =     0;
% Deleting zero error trailing at the end
err_als_vec(err_als_vec==0)=[];
%err_als_vec=[0 iter_err:iter_err:iter_err*t;err_als_vec];
%-------------------------------------------------------------------------%
end