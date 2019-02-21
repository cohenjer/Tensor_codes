% Accelerated Sparse hierarchical alternating least squares (HALS) algorithm of
% Cichocki et al. 
%
% See N. Gillis and F. Glineur, "Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization�, 
% Neural Computation 24 (4), pp. 1085-1105, 2012. 
% See http://sites.google.com/site/nicolasgillis/ 
%
% Sparse version by Jeremy E. Cohen, using a sparse nnlsHALS
% v0.1 03/10/2018
%
% [U,V,e,t] = sHALSacc(M,U,lambda,V,alpha,delta,maxiter,timelimit)
%
% Input.
%   M              : (m x n) matrix to factorize
%   lambda         : initial regularization parameter
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
% Remark. With alpha = 0, it reduces to the original HALS algorithm.  

function [U,V,e,t] = sHALS_review(M,U,lambda,V,maxiter)

% Initialization
etime = cputime; nM = norm(M,'fro')^2; 
a = 0; e = norm(M-U*V,'fro')^2; t = []; iter = 0; 

if nargin <= 4, maxiter = 100; end
if nargin <= 5, timelimit = 60; end

% Scaling, p. 72 of the thesis
A = (U'*M); B = (U'*U); j = 0;
scaling = sum(sum(V.*A))/sum(sum( B.*(V*V') )); V = V*scaling; 
% Main loop
while iter <= maxiter && e(end)>1e-6
   
    % Update of V
    if j == 1 % Do not recompute A and B at first pass
        % Use actual computational time instead of estimates rhoU
    A = (U'*M); B = (U'*U);
    end
    V = HALSupdt_sp(V,lambda,B,A,100);

    %Update of U
    A = V*M'; B = V*V'; 
    U = HALSupdt_sp(U',0,B,A,100); U = U';
    %[U,V] = col_norm(U,V);
    
    % Evaluation of the error e at time t
    if nargout >= 3
        cnT = cputime;
        %e = [e sqrt( (nM-2*sum(sum(V.*A))+ sum(sum(B.*(V*V')))) )]; 
        e = [e sqrt( (nM-2*sum(sum(A'.*U))+ sum(sum(B.*(U'*U)))) )]; 
        %e = [e norm(M - U*V,'fro')];
        etime = etime+(cputime-cnT);
        t = [t cputime-etime];
    end

    if mod(iter,100)==0
        clc
        fprintf('Sparse HALS\n')
        fprintf('Completion...%d%%\n',iter/maxiter*100)
        fprintf('Error : %d\n\n', e(end)^2)
    end
    iter = iter + 1; j = 1;
    
end

% Update of V <- HALS(M,U,V)
% i.e., optimizing min_{V >= 0} ||M-UV||_F^2 
% with an exact block-coordinate descent scheme
function V = HALSupdt(V,UtU,UtM,eit1,alpha,delta)
[r,~] = size(V); 
eit2 = cputime; % Use actual computational time instead of estimates rhoU
cnt = 1; % Enter the loop at least once
eps = 1; eps0 = 1; eit3 = 0;
while cnt == 1 || (cputime-eit2 < (eit1+eit3)*alpha && eps >= (delta)^2*eps0)
    nodelta = 0; if cnt == 1, eit3 = cputime; end
        for k = 1 : r
            deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
            V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta; 
        eit3 = cputime-eit3; 
    end
    eps = nodelta; cnt = 0; 
end

function V = HALSupdt_sp(V,lambda,UtU,UtM,itermax)
[r,~] = size(V); 
eps0 = 0;
eps = 1;
cnt = 1;
delta = 1e-2;
while cnt <= itermax && eps >= (delta)^2*eps0
    nodelta = 0; 
        for k = 1 : r
            Vkold = V(k,:);
            V(k,:) = max(((UtM(k,:)-UtU(k,:)*V)-lambda)/UtU(k,k)+V(k,:),0); % update verified
            nodelta = nodelta + norm(Vkold - V(k,:),'fro')^2; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta;
    end
    eps = nodelta; cnt = cnt+1; 
end
