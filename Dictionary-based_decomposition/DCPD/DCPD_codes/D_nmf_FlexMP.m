function [A_aux,B_aux,K_aux,err_anls_vec]=     D_nmf_FlexMP(Y,iter,B_0,K_0,D,delta,beta,opts)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,K_aux,err_anls_vec]=     D_nmf_FlexMP(Y,iter,B_0,A_0,D)
%
% Nonnegative alternating least squares algorithm for 2-way CP 
% decomposition with known dictionnary for first mode.
% Based on MP on factor A.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - A_0,B_0     : initial parameters for the factors;
% - D           : known dictionary for factor B;
% - delta       : initial value of the penalization factor;
% - beta        : minimal value of ||A-DS||/||A||;
% - opts        : options : set opts{2} to 'true' of initializing nnls with
%                 previous estimate, else the initialization is random. Set
%                 opts{3} to 'NMFinit' to refine initialization with a few iterates
%                 of HALS.
%
% Outputs:
% 
% - A_aux       : first factor (normalized)
% - K_aux       : selected atoms indexes;
% - B_aux       : second factor;
% - err_als_vec : reconstruction errors.
%
% External functions:
%
% - kr.m        : Khatri-Rao product.
%
% List of updates                 -     31/10/2016   J.E.Cohen
%                                       Creation of the file 
%                                       01/03/2017   J.E.Cohen
%                                       Introduced K, removed S
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     10;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     10^-5;%eps;
stop_err  =     eps;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
% Inner iterations for NNLS
inner     =     10;
% Dictionnary normalization for MatchP
normD     =     sqrt(sum(D.^2));
Dn = D./(repmat(normD,size(D,1),1)+1e-6); 
%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%

% Storing auxiliary variables
A_aux     =     D(:,K_0);
B_aux     =     abs(B_0);
K_aux     =     K_0;
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     Y-A_aux*B_aux';
err_anls_old    =     sqrt((res(:)')*res(:));
% Dimensions of the problem
[m,n]           =     size(Y);
% First error
err_anls_vec(1) =     err_anls_old;
fprintf('\nMP-NALS-Flexible\n')
fprintf('----\n')

ynorm   =  sqrt(sqFnorm(Y));

if strcmp(opts{3},'NMFinit')
% Few updares of the NMF algorithm HALSacc 
[A_aux,B_aux] = v2_HALSacc(Y,A_aux,B_aux',0.5,0.1,inner);
B_aux = B_aux';
K_aux             =     MatchP(A_aux,Dn);
end

disp(' Init Error ')
sqrt(sqFnorm(Y-A_aux*B_aux'))/ynorm


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
    
    % B first estimate

    if strcmp(opts{2},'true')
    B_aux    =     nnlsHALSupdt(Y,A_aux,B_aux',inner)';
    else
    B_aux    =     nnlsHALSupdt(Y,A_aux)';
    end
    
    if strcmp(opts{2},'true')
    A_aux             =     nnlsHALSupdtv2(Y',B_aux,A_aux',5*inner,D(:,K_aux)',delta)';
    else
    A_aux             =     nnlsHALSupdtv2(Y',B_aux,[],500,D(:,K_aux)',delta)';
    end
    
    % MP step
    K_aux             =     MatchP(A_aux,Dn);
    
    
    if (norm(A_aux-D(:,K_aux))/norm(A_aux) > beta) && (mod((t-1)*iter_err+i,20)==0)
        delta = delta*1.5; 
        fprintf('-----delta increased to %d-----\n',delta)
    end
    
    end
    
    % Squared norm of the error (it uses some precalculated quantities)
    err_anls          =     sqrt(sqFnorm(Y-A_aux*B_aux'))/ynorm;
    
    % Evaluation of the relative error decrease
    diff_err_filt(1:end-1)  =     diff_err_filt(2:end);
    diff_err_filt(end)      =     abs((err_anls-err_anls_old)/err_anls);
    
    % Storing the last error
    err_anls_old      =     err_anls;
    err_anls_vec(t+1) =     err_anls;
    
    if mod(t*iter_err,iter_err)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_anls,t*iter_err)
    fprintf('penalty residual : %d\n\n',norm(A_aux-D(:,K_aux))/norm(A_aux))
    end
end

    if norm(A_aux-D(:,K_aux))/norm(A_aux) > beta
        fprintf('Insufficient initial value of delta')
    end

B_aux             =     nnlsHALSupdt(Y,A_aux,B_aux',500)';
err_anls_vec(t+2) =     sqrt(sqFnorm(Y-A_aux*B_aux'))/ynorm;
fprintf('final error : %g',err_anls_vec(end))    
%-------------------------------------------------------------------------%
end