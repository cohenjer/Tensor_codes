function [A_aux,B_aux,K_aux,err_anls_vec]=     Dmult_nmf_OMP(Y,iter,B_0,A_0,D,contents,opts)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,K_aux,err_als_vec]=     D_nmf_OMP(T,iter,B_0,K_0,D,contents,opts)
%
% Nonnegative alternating least squares algorithm for 2-way CP 
% decomposition with known dictionnary for first mode Y = DS*B'.
% Based on MP on factor A=DS.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - A_0,B_0     : initial parameters for the factors;
% - D           : known dictionaries for factor B. Cell of dictionaries.
% (if possible, better to normalize each atom with l2 norm)
% - contents    : number of selected atoms in each dictionary (row vector d).
% - opts        : options :     Set opts{1} to 'NN' for nonnegative
%                constraints on all modes.
%                               Set opts{2} to 'NMFinit' to refine
%                initialization with a few iterates of HALS. 
%                               opts{3} contains the metric for the
%                assignement problem (default is 'eucl'). 
%
% Outputs:
% 
% - A_aux       : first factor (normalized)
% - K_aux       : selected atoms indexes;
% - B_aux       : second factor;
% - err_als_vec : reconstruction errors;
%
% External functions:
%
% - kr.m        : Khatri-Rao product;
% - nnlsHALSupdt: NNLS algorithm by Gillis et al (2012);
% - sqFnorm     : squared Frobenius norm for tensors;
% - MatchP      : Finds the best column sparse S in A=DS.
%
% List of updates                 -     31/10/2016   J.E.Cohen
%                                       Creation of the file 
%                                       01/03/2017   J.E.Cohen
%                                       Introduced K, removed S
%                                       04/07/2017   J.E.Cohen
%                                       multidictionary version
%-------------------------------------------------------------------------%

% Setting the options if not provided
if nargin<7
    opts{1} = 'NN'; 
    opts{2} = 'NMFinit';
    opts{3} = 'eucl';
end

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     10;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     10^-5;
stop_err  =     eps;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
% Inner iterations for NNLS
inner     =     10;
%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%

% Storing auxiliary variables
A_aux     =     abs(A_0);
B_aux     =     abs(B_0);
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     Y-A_aux*B_aux';
err_anls_old    =     sqrt((res(:)')*res(:));
% First error
err_anls_vec(1) =     err_anls_old;
fprintf('\n M2PALS Running\n')
fprintf('----\n')

ynorm=sqrt(sqFnorm(Y));

if strcmp(opts{2},'NMFinit')
% Few updares of the NMF algorithm HALSacc 
[A_aux,B_aux] = v2_HALSacc(Y,A_aux,B_aux',0.5,0.1,inner);
B_aux = B_aux';
end


% disp(' Init Error ')
% sqrt(sqFnorm(Y-D(:,MatchP(A_aux,Dn))*B_aux'))/ynorm
% sqrt(sqFnorm(Y-A_aux*B_aux'))/ynorm


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

    % MP step
    [alloc,indexes,A_aux]    =     MatchP_multiD(A_aux,D,contents,opts{3});
    K_aux             =     [alloc;indexes];
    A_old             =     A_aux; % For cost function
    
    % B first estimate
    if strcmp(opts{1},'NN')
        B_aux    =     nnlsHALSupdt(Y,A_aux,B_aux',inner)';
    else
        B_aux    =     (A_aux\Y)';
    end
    
    % A update
    
    if strcmp(opts{1},'NN')
    A_aux             =     nnlsHALSupdt(Y',B_aux,A_aux',inner)';
    else
    A_aux             =     Y/(B_aux');
    end
    

    end
    
    % Squared norm of the error (it uses some precalculated quantities)
    err_anls          =     sqrt(sqFnorm(Y-A_old*B_aux'))/ynorm;
    A_aux             =     A_old; % modification from original output which was A_aux
    
    % Evaluation of the relative error decrease
    diff_err_filt(1:end-1)  =     diff_err_filt(2:end);
    diff_err_filt(end)      =     abs((err_anls-err_anls_old)/err_anls);
    
    % Storing the last error
    err_anls_old      =     err_anls;
    err_anls_vec(t+1) =     err_anls;
    
    if mod(t,5)  ==    0
    % % Printing the iterate error at each 50 iterates
    fprintf('error : %g  \tit : %d\n',err_anls,t*iter_err)
    end
end

if strcmp(opts{1},'NN')
B_aux             =     nnlsHALSupdt(Y,A_aux,B_aux',500)';
end
err_anls_vec(t+2)          =     sqrt(sqFnorm(Y-A_aux*B_aux'))/ynorm;
% fprintf('final error : %g\n',err_anls_vec(end))
%-------------------------------------------------------------------------%
end