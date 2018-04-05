function V = nnlsHALSupdtv2(M,U,V,maxiter,Vtarget,lambda)

% Computes an approximate solution of the following nonnegative least 
% squares problem (NNLS)  
%
%           min_{V >= 0} ||M-UV||_F^2 + lambda ||V - Vtarget||_F^2   (1)
% 
% with an exact block-coordinate descent scheme. 
%
% See N. Gillis and F. Glineur, Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization, 
% Neural Computation 24 (4): 1085-1105, 2012.
% 
%
% ****** Input ******
%   M  : m-by-n matrix 
%   U  : m-by-r matrix
%   V  : r-by-n initialization matrix 
%        default: one non-zero entry per column corresponding to the 
%        clostest column of U of the corresponding column of M 
%   maxiter: upper bound on the number of iterations (default=500).
%
%   *Remark. M, U and V are not required to be nonnegative. 
%
% ****** Output ******
%   V  : an r-by-n nonnegative matrix that solves (1) approximately

if nargin <= 3
    maxiter = 500;
end
[m,n] = size(M); 
[m,r] = size(U); 
UtU = U'*U; 
UtM = U'*M;

if nargin <= 4
    Vtarget = zeros(r,n);
    lambda = 0; 
end

if nargin <= 2 || isempty(V) 
    V = U\M; % Least Squares
    V = max(V,0); 
    % Scaling 
    alpha = sum(sum( (U'*M).*V ))./sum(sum( (U'*U).*(V*V'))); 
    V = alpha*V; 
end

delta = 1e-6; % Stopping condition depending on evolution of the iterate V: 
              % Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F 
              % where V^{k} is the kth iterate. 
eps0 = 0; cnt = 1; eps = 1; 
while eps >= delta*eps0 && cnt <= maxiter %Maximum number of iterations
    nodelta = 0; if cnt == 1, eit3 = cputime; end
    Vold = V; 
        for k = 1 : r
            V(k,:) = max( 0 , ( UtM(k,:)-UtU(k,:)*V+UtU(k,k)*V(k,:) + lambda*Vtarget(k,:) ) / (UtU(k,k)+lambda)); 
            if V(k,:) == 0, 
                V(k,:) = 1e-16*max(V(:));  % safety procedure
            end
        end
    nodelta = norm(V-Vold,'fro'); 
    if cnt == 1
        eps0 = nodelta; 
    end
    eps = nodelta; 
    cnt = cnt + 1; 
end
end % of function nnlsHALSupdt