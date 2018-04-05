% Accelerated hierarchical alternating least squares (HALS) algorithm of
% Cichocki et al. 
%
% See N. Gillis and F. Glineur, "Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization”, 
% Neural Computation 24 (4), pp. 1085-1105, 2012. 
% See http://sites.google.com/site/nicolasgillis/ 
%
% [U,V,e,t] = v2_HALSacc(M,U,V,alpha,delta,maxiter,timelimit)
%
% Input.
%   M              : (m x n) matrix to factorize
%   (U,V)          : initial matrices of dimensions (m x r) and (r x n)
%   alpha          : nonnegative parameter of the accelerated method
%                    (alpha=0.5 seems to work well)
%   delta          : parameter to stop inner iterations when they become
%                    inneffective (delta=0.1 seems to work well). 
%   maxiter        : maximum number of iterations
%   timelimit      : maximum time alloted to the algorithm
%
% Output.
%   (U,V)    : nonnegative matrices s.t. UV approximate M
%   (e,t)    : error and time after each iteration, 
%               can be displayed with plot(t,e)
%
% Remark 1. With alpha = 0, it reduces to the original HALS algorithm.  
% Remark 2. Compared to the code HALSacc from the original paper above, 
% v2_HALSacc updates U and V a number of times that depends on the
% theoretical computational cost, as opposed to the CPU time for HALSacc. 
% The advantage is that it makes the algorithm deterministic hence the 
% numerical experiments are reproducible. 

function [U,V,e,t] = v2_HALSacc(M,U,V,alpha,delta,maxiter,timelimit)

% Initialization
etime = cputime; nM = norm(M,'fro')^2; 
[m,n] = size(M); [m,r] = size(U);
if issparse(M), 
    K = sum(M(:) ~= 0); 
else K = m*n; 
end
rhoU = 1+(K+n*r)/(m*(r+1)); 
rhoV = 1+(K+m*r)/(n*(r+1)); 
maxiterU = floor(1+rhoU*alpha); 
maxiterV = floor(1+rhoV*alpha); 

a = 0; e = []; t = []; iter = 0; 

if nargin <= 3, alpha = 0.5; end
if nargin <= 4, delta = 0.1; end
if nargin <= 5, maxiter = 100; end
if nargin <= 6, timelimit = 60; end

% Scaling, p. 72 of the thesis
A = M*V'; B = V*V'; 
j = 0;
scaling = sum(sum(A.*U))/sum(sum( B.*(U'*U) )); 
U = U*scaling; 
% Main loop
while iter <= maxiter && cputime-etime <= timelimit
    % Update of U
    if j == 1, % Do not recompute A and B at first pass
        % Use actual computational time instead of estimates rhoU
        eit1 = cputime; A = M*V'; B = V*V'; eit1 = cputime-eit1; 
    end
    j = 1; 
    U = HALSupdt(U',B',A',maxiterU,alpha,delta); U = U';
    
    % Update of V
    A = (U'*M); B = (U'*U); 
    V = HALSupdt(V,B,A,maxiterV,alpha,delta); 
    % Evaluation of the error e at time t
    if nargout >= 3
        cnT = cputime;
        e = [e sqrt( (nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) )]; 
        etime = etime+(cputime-cnT);
        t = [t cputime-etime];
    end
    iter = iter + 1; j = 1; 
end

% Update of V <- HALS(M,U,V)
% i.e., optimizing min_{V >= 0} ||M-UV||_F^2 
% with an exact block-coordinate descent scheme
function V = HALSupdt(V,UtU,UtM,maxiter,alpha,delta)

[r,n] = size(V); 
eps = 1; eps0 = 1; 
j = 1; 
while j <= maxiter && eps >= (delta)^2*eps0
    nodelta = 0;
    for k = 1 : r
        deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
        V(k,:) = V(k,:) + deltaV;
        nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
        if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
    end
    if j == 1
        eps0 = nodelta; 
    end
    eps = nodelta; 
    j = j+1; 
end