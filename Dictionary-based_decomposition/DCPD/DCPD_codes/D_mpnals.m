function [A_aux,B_aux,C_aux,K_aux,err_anls_vec]=     D_mpnals(Y,iter,A_0,B_0,C_0,D,opts)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,K_aux,err_anls_vec]=     D_mpnals(Y,iter,A_0,B_0,C_0)
%
% Alternating least squares algorithm for 3-way CP 
% decomposition with known dictionnary for second mode. May use
% non-negative least squares for non-negative data.
% Factors A and B are normalized with the L_2 norm at each iterate.
% Based on OMP on factor C. Using acc-HALS from Gillis et. al. for the
% nonnegative least squares update.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - A_0,C_0     : initial parameters for the factors;
% - S_0         : initial score factor;
% - D           : known dictionary for factor B.
% - opts        : options. 
%       opts{1} = 'NN' for non-negative updates
%       opts{2} = 'true' for initializing nnls with previous estimate (NN case only).
%
% Outputs:
% 
% - A_aux       : first factor;
% - S_aux       : second factor scores;
% - C_aux       : third factor;
% - err_als_vec : reconstruction errors 
%
% External functions:
%
% - kr.m        : Khatri-Rao product.
% - nnlsHALSupdt: Non-negative least squares algorithm.
% - MatchP      : finds the best atoms in the dictionnary
%
% List of updates                 -     23/02/2016   J.E.Cohen
%                                       Creation of the file 
%                                 -     08/02/2013   J.E.Cohen
%                                       Added non-negativity and
%                                       initialization options.
%                                 -     01/03/2017   J.E.Cohen
%                                       removed S, introduced K
%                                 -     04/07/2017   J.E.Cohen
%                                       dictionary normalization
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     100;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     10^-8;
stop_err  =     eps;
% Size of the tensor
dim       =     size(Y);
% Index for n-mode matricization
J         =     prod(dim)./dim;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
% Inner iterations for NNLS
inner     =     20;
%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%
% Three different unfoldings of the tensor
Y_mode    =     cell(3);
% Data tensor unfolding
Y_mode{1} =     reshape(permute(Y,[1,2,3]),dim(1),J(1));  
Y_mode{2} =     reshape(permute(Y,[2,1,3]),dim(2),J(2));
Y_mode{3} =     reshape(permute(Y,[3,1,2]),dim(3),J(3));
% Storing auxiliary variables
A_aux     =     A_0;
B_aux     =     B_0;
C_aux     =     C_0;
C_aux     =     C_aux.*repmat(sqrt(sum(B_aux.^2)),size(C_aux,1),1);
B_aux     =     col_norm(B_aux);
% Dictionnary normalization for MatchP
normD     =     sqrt(sum(D.^2));
Dn = D./(repmat(normD,size(D,1),1)+1e-6); 
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     Y_mode{1}-(A_aux)*transpose(kr(C_aux,B_aux));
err_anls_old    =     (res(:)')*res(:);
% Squared norm of Y
snorm_Y         =     Y_mode{1}(:)'*Y_mode{1}(:);
% Vector of squared norms of the error
%err_anls_vec    =     zeros(1,iter+1);
% First error
err_anls_vec(1) =     err_anls_old;
fprintf('\nMPALS\n')
fprintf('----\n')

%-------------------------------------------------------------------------%

%--------------------------Algorithm iterate------------------------------%
t         =     0;
while     (mean(diff_err_filt)>stop_diff &&  t<iter ...
          &&  err_anls_old>stop_err   )
    % Iteration counter (counts/iter_err)
    t=    t+1;
    
    % ANLS (iter_err iterates)
    for   i           =     1:iter_err
    
    % Least squares solution for each factor
    
    % A update
    
if strcmp(opts{1},'NN')
    if strcmp(opts{2},'true')
        A_aux             =     nnlsHALSupdt(Y_mode{1}',kr(C_aux,B_aux),A_aux',inner)';
    else
        A_aux             =     nnlsHALSupdt(Y_mode{1}',kr(C_aux,B_aux))';
    end
else
        Y_BC              =     Y_mode{1}*kr(C_aux,B_aux);
        temp_BC           =     (C_aux'*C_aux).*(B_aux'*B_aux);
        A_aux             =     Y_BC/temp_BC;
end
    A_aux             =     col_norm(A_aux);
    
    % B first estimate
if strcmp(opts{1},'NN')    
    if strcmp(opts{2},'true')
    B_aux             =     nnlsHALSupdt(Y_mode{2}',kr(C_aux,A_aux),B_aux',inner)';
    else
    B_aux             =     nnlsHALSupdt(Y_mode{2}',kr(C_aux,A_aux))';
    end
else
    Y_AC              =     Y_mode{2}*kr(C_aux,A_aux);
    temp_AC           =     (C_aux'*C_aux).*(A_aux'*A_aux);
    B_aux             =     Y_AC/temp_AC;
end
    norm_B            =     sum(B_aux.^2);
    B_aux             =     B_aux.*repmat(1./sqrt(norm_B),dim(2),1);
    % OMP step
    K_aux             =     MatchP(B_aux,Dn);
    % B final estimate
    B_aux             =     D(:,K_aux);

    % C update
if strcmp(opts{1},'NN')
    if strcmp(opts{2},'true')
    C_aux             =     nnlsHALSupdt(Y_mode{3}',kr(B_aux,A_aux),C_aux',inner)';
    else
    C_aux             =     nnlsHALSupdt(Y_mode{3}',kr(B_aux,A_aux))';
    end
else
    Y_AB              =     Y_mode{3}*kr(B_aux,A_aux);
    temp_AB           =     (B_aux'*B_aux).*(A_aux'*A_aux);
    C_aux             =     Y_AB/temp_AB;
end
   
    end
    
    % Squared norm of the error (it uses some precalculated quantities)
    temp_BA           =     (B_aux'*B_aux).*(A_aux'*A_aux);
    Y_BA              =     Y_mode{3}*kr(B_aux,A_aux);
    err_anls          =     snorm_Y-2*sum(sum(Y_BA.*C_aux))+...
                            sum(sum((C_aux'*C_aux).*temp_BA));
    
    % Evaluation of the relative error decrease
    diff_err_filt(1:end-1)  =     diff_err_filt(2:end);
    diff_err_filt(end)      =     abs((err_anls-err_anls_old)/err_anls);
    
    % Storing the last error
    err_anls_old      =     err_anls;
    err_anls_vec(t+1) =     err_anls;
    
    if mod(t*iter_err,iter_err)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_anls,t*iter_err)
    end
end

%-------------------------------------------------------------------------%
end
