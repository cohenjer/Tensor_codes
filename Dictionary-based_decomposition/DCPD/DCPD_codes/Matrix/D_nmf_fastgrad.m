function [A_aux,B_aux,K_aux,err_anls_vec]=     D_nmf_fastgrad(Y,iter,A_0,B_0,D,opts)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,err_als_vec]=      D_nmf_fastgrad(Y,iter,A_0,S_0,D)
%
% Nesterov's Fast Gradient for NMF with known dictionnary for second mode. The projection can
% be either alternated for l1 and l2 norms or simultaneous.
% Factors A is normalized with the L_2 norm at each iterate.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - A_0         : initial parameters for the factors;
% - B_0         : initial spectra;
% - D           : known dictionary for factor B (normalized).
% - opts        : options : set opts{2} to 'true' of initializing nnls with
%                 previous estimate, else the initialization is random.
%
% Outputs:
% 
% - A_aux       : first factor;
% - B_aux       : second factor;
% - K_aux       : selected atoms indexes;
% - err_als_vec : reconstruction errors
%
% External functions:
%
% - kr.m        : Khatri-Rao product;
% - col_norm    : l2 normalisation of columns;
% - simplex_proj: projection on the probability simplex.
%
% List of updates                 -     04/01/2017   J.E.Cohen
%                                       Creation of the file 
%                                       01/03/2017   J.E.Cohen
%                                       Introduced K
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     1;
% Number of inner iterations
iter_inner = 10; % Note : reduce if algorithm is diverging.
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Dimensions of the problem
[m,n]           =     size(Y);
% Stopping criteria parameters
stop_diff =     10^-5;%10^-4
stop_err  =     eps;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
% NNLS iteration number
inner     =     30;
% l1 regularization parameter
lambda_max    =     10;%0.05

%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%
% Storing auxiliary variables
A_aux     =     abs(A_0);
B_aux     =     abs(B_0);
K_aux     =     MatchP(B_aux,D);
S_aux     =     0.2*rand(size(D,2),size(B_aux,2));
for j=1:size(B_aux,2)
S_aux(K_aux(j),j)     =     1;
end
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     Y-A_aux*B_aux';
err_anls_old    =     sqrt((res(:)')*res(:));
% Squared norm of Y
norm_Y         =     sqrt(Y(:)'*Y(:));
% Vector of squared norms of the error
%err_anls_vec    =     zeros(1,iter+1);
% First error
err_anls_vec(1) =     err_anls_old/norm_Y;
fprintf('\n Alternating Fast Gradient \n')
fprintf('----\n')

if strcmp(opts{3},'NMFinit')
% Few updares of the NMF algorithm HALSacc 
[A_aux,B_aux] = v2_HALSacc(Y,A_aux,B_aux',0.5,0.1,inner);
B_aux = B_aux';
B_aux    =     B_aux.*repmat(sqrt(sum(A_aux.^2)),n,1);
A_aux    =     A_aux.*repmat(1./sqrt(sum(A_aux.^2)),m,1);
end


disp(' Init Error ')
sqrt(sqFnorm(Y-A_aux*D(:,MatchP(B_aux,D))'))/norm_Y


% Accelerated computations

%Gram = D'*D;
[~,stepD,~]=svds(D,1);stepD = stepD^2;
%DT   = D'*Y';

% Others
alpha_S_new = 0.5;%0.5

%-------------------------------------------------------------------------%

%--------------------------Algorithm iterate------------------------------%
t         =     0;

while     (mean(diff_err_filt)>stop_diff &&  t<iter ...
          &&  err_anls_old>stop_err   )
    % Iteration counter (counts/iter_err)
    t=    t+1;
    
    %if nnz(S_aux)>R
    lambda = lambda_max*((t-1)/iter);
    %end
    % Gradient iterates
    for   i           =     1:iter_err
    
        
    % S gradient        
    
        temp_A_aux    =  A_aux'*A_aux;
        [~,stepB,~]   =  svds(temp_A_aux,1);
        stepS         =  stepD*stepB;
        temp_DY       =  D'*(Y'*A_aux);
    
        for j=1:iter_inner
            
        S_old    =  S_aux;
        GradS    = -temp_DY+D'*(D*(S_aux*(temp_A_aux)))+lambda*ones(size(S_aux));
    
            % Direction and acceleration
            
        S_aux    =  max(S_aux - 1/stepS*GradS,0);
        S_aux    =  col_norm(S_aux);
        S_aux(S_aux<10^(-16))=0;
 
        
        alpha_S_old = alpha_S_new;
        alpha_S_new = 1/2*(-alpha_S_old^2+sqrt(alpha_S_old^4+4*alpha_S_old^2));
        beta_S      = alpha_S_old*(1-alpha_S_old)/(alpha_S_old^2+alpha_S_new);
        S_aux       = S_aux + beta_S*(S_aux-S_old); 
        
        % Visualization
        pause(0.05)
        spy(S_aux(1:100,:))
        
        end
                  
        B_aux  =  D*S_aux;
        
            % A update    
       
        if strcmp(opts{2},'true')
        A_aux             =     nnlsHALSupdt(Y',B_aux,A_aux',inner)';
        else
        A_aux             =     nnlsHALSupdt(Y',B_aux)';        
        end
        
    end
    
    % Squared norm of the error (of the previous iterate)
    err_anls          =     sqrt(sqFnorm(Y-A_aux*D(:,MatchP(B_aux,D))'))/norm_Y;
    
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

    K_aux     = MatchP(B_aux,D);
    B_aux     = D(:,K_aux); 
    
 % To be used for published version, or run in independent code   
%A_aux             =     nnlsHALSupdt(Y',D(:,K_aux),A_aux',500)';
%err_anls_vec(t+2)          =     sqrt(sqFnorm(Y-A_aux*D(:,K_aux)'))/norm_Y;
%fprintf('final error : %g',err_anls_vec(end))
    
%-------------------------------------------------------------------------%
end