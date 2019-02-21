function [A,B,err] = k_s_nmf(M,k,R,iter,A0,B0)
% [A,B,B_aux,err] = k_s_nmf(M,k,R,iter,A0,B0)
% This function computes the k-sparse nonnegative matrix factorization of
% data matrix M:
%               M = AB', ||b_i||_0 < k+1
% K-sparisty is imposed on the columns of matrix B, so that
% at most k components are active at each sample along the second mode. 
%
% Implemented with brute force search for enforcing k-sparsity.
%
% Inputs :
% - M : n x m nonnegative data matrix
% - k : number of nonzero coefficients in columns of B
% - R : rank of the factorization 
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
% 16/01/2018 J. E.Cohen
% Creation of the file
% 17/01/2018 J. E.Cohen
% Vectorizing the B computation, changing output visuals
% 24/01/2018 J. E.Cohen
% Added stopping criterion

A = A0;
B = B0;
normM = norm(M,'fro');
err(1) = norm(M-A*B,'fro')/normM;  err(2) = inf;
tol = 10^-5; % minimal relative reconstruction error
t = 0; % loop counter
iter_nn = 50; % Number of HALS iterations for nnls
iter_perm = nchoosek(R,k); % Number of permutations for brute force
perms_B = make_perm(R,k,1); % Obtain set of permuted indexes (Homemade)

% Main loop
while t<iter && abs((err(end)-err(end-1))/err(end-1))>tol && err(end) > 10^-30
t = t+1;
clc
fprintf('K-sparse NMF running...\n')
fprintf('Completion... ')
fprintf('%d%%  \n',round(100*t/iter))
if t==1
fprintf('Error : %d  ',err(1))
else
fprintf('Error : %d  ',err(end))    
end

% B estimate with brute force search
err_i = zeros(iter_perm,size(M,2));
Bi = cell(iter_perm,1);
for j=1:iter_perm % loop over chosen k coefficients
   indexes = perms_B(j,:);
   % Computing LS with k coefficients
   Bi{j} = nnlsHALSupdt(M,A(:,indexes),B(indexes,:),200);%iter_nn);
   % Computing error (row vector)
   err_i(j,:) = sum((M-A(:,indexes)*Bi{j}).^2);
end
% Finding best permutations for all coefficient vector
[~,best_set] = min(err_i);
% Affecting best results to B
for i=1:size(M,2)
    B(:,i) = 0;
    B(perms_B(best_set(i),:),i) = Bi{best_set(i)}(:,i);
end

% A estimate using nonnegative least squares
A = nnlsHALSupdt(M',B',A',iter_nn)';
%A = A.*repmat(1./sqrt(sum(A.^2)),size(A,1),1);
[A,B] = col_norm(A,B);

% error computation
err(t+1) = norm(M-A*B,'fro')/normM;


end


clc
fprintf('K-sparse NMF running...\n')
fprintf('Completion... ')
fprintf('100%%  ')

function [ M,N ] = col_norm( M,N )

[~,n] = size(M);

    if nargin==1
norms = sqrt(sum(M.^2));
for r=1:n
    if norms(r)>0
        M(:,r) = M(:,r)/norms(r);
    end
end
N=[];
    else
norms = sqrt(sum(M.^2));
for r=1:n
    if norms(r)>0
        M(:,r) = M(:,r)/norms(r);
        N(:,r) = N(:,r)*norms(r);
    end
end

    end
end




end

