function [A_aux,B_aux,C_aux,err_anls_vec,temps]=     anls(Y,iter,A_0,B_0,C_0)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,err_anls_vec]=     anls(Y,iter,A_0,B_0,C_0)
%
% Alternating nonnegative least squares algorithm for 3-way CP
% decomposition. Factors A and B are normalized with the L_1 norm at each
% iterate.
%
% Inputs:
% 
% - Y           : data block;
% - iter        : maximum number of iterates;
% - R           : number of components;
% - A_0,B_0,C_0 : initial parameters for the factors.
%
% Outputs:
% 
% - A_aux       : first factor;
% - B_aux       : second factor;
% - C_aux       : third factor;
% - err_anls_vec: iterations (Row 1) x reconstruction errors (Row 2).
%
% External functions:
%
% - kr.m              : Khatri-Rao product.
%
% List of updates                 -     03/01/2015  -     R. C. Farias
%                                       Creation of the file
%                                 -     20/07/2015  -     R. C. Farias
%                                       Update with fast right inverse
%-------------------------------------------------------------------------%

%-------------------------------Parameters--------------------------------%
% Number of iterations for each evaluation of the reconstruction error
iter_err  =     1;
% Size of the filter for the relative error decrease
n_diff_err=     3;
% Stopping criteria parameters
stop_diff =     eps;
stop_err  =     eps;
% Size of the tensor
dim       =     size(Y);
% Index for n-mode matricization
J         =     prod(dim)./dim;
% Maximum number of iterations (/iter_err)
iter      =     ceil(iter/iter_err);
%-------------------------------------------------------------------------%

%-------------Missing values ---------------
% Setting unknown values to 0
vecy      =     Y(:);
M         =     isnan(vecy);
vecy(M)   =     min(vecy);
Y         =     reshape(vecy,dim);


%-----------------------Algorithm initialization--------------------------%
% Three different unfoldings of the tensor
Y_mode    =     cell(3);
% Data tensor unfolding
Y_mode{1} =     reshape(permute(Y,[1,2,3]),dim(1),J(1));  
Y_mode{2} =     reshape(permute(Y,[2,1,3]),dim(2),J(2));
Y_mode{3} =     reshape(permute(Y,[3,1,2]),dim(3),J(3));
% Storing auxiliary variables
A_aux     =     abs(A_0);
B_aux     =     abs(B_0);
C_aux     =     abs(C_0);
% Filter for the relative error decrease
diff_err_filt   =     1000*ones(1,n_diff_err);
% Initial squared norm of the error
res             =     Y_mode{1}-(A_aux)*transpose(kr(C_aux,B_aux));
err_anls_old    =     (res(:)')*res(:);
% Squared norm of Y
snorm_Y         =     Y_mode{1}(:)'*Y_mode{1}(:);
% Vector of squared norms of the error
err_anls_vec    =     zeros(1,iter+1);
% First error
err_anls_vec(1) =     err_anls_old;
fprintf('\nANLS\n')
fprintf('----\n')
%-------------------------------------------------------------------------%

%--------------------------Algorithm iterate------------------------------%
t         =     0;
temps(1)  =     1;

while     (mean(diff_err_filt)>stop_diff &&  t<iter ...
          &&  err_anls_old>stop_err   )
    % Iteration counter (counts/iter_err)
    t=    t+1;
    ticker=tic;
    % ANLS (iter_err iterates)
    for   i           =     1:iter_err
    
    % Least squares solution for each factor
    
    % A update
    temp_CB           =     (C_aux'*C_aux).*(B_aux'*B_aux);
    A_aux             =     (Y_mode{1}*kr(C_aux,B_aux))/temp_CB;
    A_aux(A_aux<0)    =     0;
    norm_A            =     sqrt(sum(A_aux.^2));
    A_aux             =     A_aux.*repmat(1./norm_A,dim(1),1);
    
    % B update
    temp_CA           =     (C_aux'*C_aux).*(A_aux'*A_aux);
    B_aux             =     (Y_mode{2}*kr(C_aux,A_aux))/temp_CA;
    B_aux(B_aux<0)    =     0;
    B_aux             =     B_aux.*repmat(1./sqrt(sum(B_aux.^2)),dim(2),1);

    % C update
    temp_BA           =     (B_aux'*B_aux).*(A_aux'*A_aux);
    Y_BA              =     Y_mode{3}*kr(B_aux,A_aux);
    C_aux             =     Y_BA/temp_BA;
    C_aux(C_aux<0)    =     0;
    end
temps(t+1)    =    toc(ticker);       
    % Squared norm of the error (it uses some precalculated quantities)
    err_anls          =     Y(:)'*Y(:)-2*sum(sum(Y_BA.*C_aux))+...
                            sum(sum((C_aux'*C_aux).*temp_BA));
    
    % Evaluation of the relative error decrease
    diff_err_filt(1:end-1)  =     diff_err_filt(2:end);
    diff_err_filt(end)      =     abs((err_anls-err_anls_old)/err_anls);
    
    % Storing the last error
    err_anls_old      =     err_anls;
    err_anls_vec(t+1) =     err_anls;
    
    if mod(t*iter_err,100)  ==    0
    % Printing the iterate error at each 100 iterates
    fprintf('error : %g  \tit : %d\n',err_anls,t*iter_err)
    end
 
    Y_est     =     MultProd(eyeT(size(A_0,2)),A_aux,B_aux,C_aux);
       
    % Estimation of the missing data

    vecy_est  =     Y_est(:);
    vecy(M)   =     vecy_est(M);
    Y         =     reshape(vecy,dim);
    
    
end
A_aux(isnan(A_aux))=1;
B_aux(isnan(B_aux))=1;
C_aux(isnan(C_aux))=1;
% Deleting zero error trailing at the end
err_anls_vec(err_anls_vec==0)     =     [];
%-------------------------------------------------------------------------%
end