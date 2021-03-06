function [A,B,err] = as_k_s_nmf(M,k,iter,A0,B0)
% [A,B,B_aux,err] = as_k_s_nmf(M,k,iter,A0,B0)
% This function computes the k-sparse nonnegative matrix factorization of
% data matrix M:
%               M = AB', ||b_i||_0 < k+1
% K-sparisty is imposed on the columns of matrix B, so that
% at most k components are active at each sample along the second mode. 
%
%
% Implemented with early stopped active set (sNNLS in [ref])
%
% Inputs :
% - M : n x m nonnegative data matrix
% - k : number of nonzero coefficients in columns of B
% - iter : maximal number of iterations (optional)
% - A0 : initial factor matrix A (random if not provided)
% - B0 : initial factor matrix B (random if not provided)
%
% Outputs :
% - A : estimated factor A
% - B : estimated k-sparse factor B
% - err : reconstruction error (relative)
%
% Updates :
% 03/10/2018 J. E.Cohen    --   Creation of the file
% 

A = A0;
B = B0;
normM = norm(M,'fro');
err(1) = norm(M-A*B,'fro')/normM;  err(2) = inf;
tol = 10^-5; % minimal relative reconstruction error
t = 0; % loop counter
iter_nn = 50; % Number of HALS iterations for nnls
options.MAX_ITER = k+1; % This ensures at most k coefficients in columns of B

% Main loop
while t<iter && abs((err(end)-err(end-1))/err(end-1))>tol && err(end) > 10^-15
t = t+1;
clc
fprintf('Active-set K-sparse NMF running...\n')
fprintf('Completion... ')
fprintf('%d%%  \n',round(100*t/iter))
if t==1
fprintf('Error : %d  ',err(1))
else
fprintf('Error : %d  ',err(end))    
end

% B estimate with early stop active set
B = nnlsm_activeset( A, M, options); 

% A estimate using nonnegative least squares
A = nnlsHALSupdt(M',B',A',iter_nn)';
%A = A.*repmat(1./sqrt(sum(A.^2)),size(A,1),1);
[A,B] = col_norm(A,B);

% error computation
err(t+1) = norm(M-A*B,'fro')/normM;


end


clc
fprintf('active-set NNOMP K-sparse NMF running...\n')
fprintf('Completion... ')
fprintf('100%%  ')

end

