function [V,e] = FGMfcnls(M,U,V,maxiter) 

% Fast gradient method to solve least squares on the unit simplex.  
% See Nesterov, Introductory Lectures on Convex Optimization: A Basic 
% Course, Kluwer Academic Publisher, 2004. 
% 
% This code solves: 
% 
%             min_{V(:,j) in Delta, forall j}  ||M-UV||_F^2, 
%
% where Delta = { x | sum x_i <= 1, x_i >= 0 for all i }.
%  
% See also Appendix A in N. Gillis, Successive Nonnegative Projection 
% Algorithm for Robust Nonnegative Blind Source Separation, arXiv, 2013. 
% 
%
% [V,e] = FGMfcnls(M,U,V,maxiter) 
%
% ****** Input ******
% M      : m-by-n data matrix
% U      : m-by-r basis matrix
% V      : initialization for the fast gradient method 
%          (optional, use [] if none)
% maxiter: maximum numbre of iterations (default = 500). 
%
% ****** Output ******
% V      : V(:,j) = argmin_{x in Delta}  ||M-Ux||_F^2 forall j. 
% e      : e(i) = error at the ith iteration

[m,n] = size(M); 
[m,r] = size(U); 
% Initialization of V 
if nargin <= 2 || isempty(V)
    V = zeros(r,n); 
    for i = 1 : n
        % Distance between ith column of M and columns of U
        disti = sum( (U - repmat(M(:,i),1,r)).^2 ); 
        [a,b] = min(disti); 
        V(b,i) = 1; 
    end
end
if nargin <= 3
    maxiter = 500; 
end
nM = norm(M,'fro')^2;
% Hessian and Lipschitz constant 
UtU = U'*U; 
L = norm(UtU,2); 
% Linear term 
UtM = U'*M; 
nM = norm(M,'fro')^2; 
alpha0 = 0.05; % Parameter, can be tuned. 
alpha(1) = alpha0;
V = SimplexProj( V ); % Project initialization onto the simplex
Y = V; % second sequence
i = 1; 
% Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F
delta = 1e-6;
eps0 = 0; eps = 1;  
while i <= maxiter && eps >= delta*eps0
    % Previous iterate
    Vp = V; 
    % FGM Coefficients  
    alpha(i+1) = ( sqrt(alpha(i)^4 + 4*alpha(i)^2 ) - alpha(i)^2) / (2); 
    beta(i) = alpha(i)*(1-alpha(i))/(alpha(i)^2+alpha(i+1)); 
    % Projected gradient step from Y
    V = SimplexProj( Y - (UtU*Y-UtM) / L );
    % `Optimal' linear combination of iterates
    Y = V + beta(i)*(V-Vp); 
    % Error 
    e(i) = nM - 2*sum(sum(V.*(UtM))) + sum(sum((UtU).*(V*V'))); 
    % Restart: fast gradient methods do not guarantee the objective
    % function to decrease, a good heursitic seems to restart whenever it
    % increases although the global convergence rate is lost! This could
    % be commented out. 
    if i >= 2 && e(i) > e(i-1)
        Y = V; 
    end 
    if i == 1
        eps0 = norm(V-Vp,'fro'); 
    end
    eps = norm(V-Vp,'fro'); 
    i = i + 1; 
end  