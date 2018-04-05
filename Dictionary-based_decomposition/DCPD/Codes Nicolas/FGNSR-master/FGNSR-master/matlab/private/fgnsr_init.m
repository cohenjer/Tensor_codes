function [X, p, mu] = fgnsr_init(M, r, p, init_K, init_X)
% [X, p, mu] = fgnsr_init(M, r, p, init_K, init_X)
%
% Outputs:
%   X    -- initial value for X
%   p    -- scaled objective for the trace
%   mu   -- scaling factor used for p
%
%   init_K  --  method for initial K, affects the scaling of p
%               0: random choice
%               1: use FastSepNMF (default)
%               2: use KKT (and FastSepNMF)
%
%   init_X  --  0: set X to zero
%               1: set X according to optimal H for the choice of K (def.)
%               2: set X to r/n * eye(n)
%               3: set X to the identity

if nargin < 5 || isempty(init_X)
    init_X = 1;
end

if nargin < 4 || isempty(init_K)
    init_K = 1;
end

n = size(M, 2);

if init_K == 0
    K = randperm(n);
    K = K(1:r);
elseif init_K == 1 || init_K == 2
    K = sort(FastSepNMF(M,r,0));
else
    assert(false);
end

H = nnlsHALSupdt(M, M(:,K)); 
X = zeros(n,n); 
X(K,:) = H;

% Evaluate approximation error for this tentative X
nM = norm(M,'fro')^2; 
U = M(:,K); 
% err equals to norm(M - M*X, 'fro')^2
err = (nM-2*sum(sum(H.*(U'*M))) + sum(sum((U'*U).*(H*H')))); 

% If X is a very good solution, we still don't want the trace objective
% to be scaled near to the zero vector.  We heuristically add some epsilon.
e = err + 1e-3;

if init_K <= 1
    %p'*diag(X)
    ep = sum(p(K)); % = p'*diag(X)
    mu = e/ep;
    p = p*e/ep; % so that ||M-MX||_F^2 = p'*diag(X)
elseif init_K == 2
    a = diag(-M'*M + M'*M*X);
    muit = linspace(0,max(abs(a)),n);
    Kbar = setdiff(1:n,K);
    evalf = norm( a(K) )^2 + norm( max(-a(Kbar),0) )^2;
    evalfcprec = +inf;
    i = 1;
    while evalfcprec > evalf
        mu = muit(i);
        evalfcprec = evalf;
        evalf = norm( a(K) + muit(i+1) )^2 + norm( max(-a(Kbar) - muit(i+1),0) )^2;
        i = i+1; 
    end
    p = p*mu; 
end



if init_X == 0
    % reset X to zero
    X(:) = 0.0;
elseif init_X == 1
    % do nothing, X is initialized with H
elseif init_X == 2
    % equal weights on the diagonal having optimal trace
    X = diag(r/n * ones(n,1));
elseif init_X == 3
    % start with the identity
    X = eye(n);
else
    % Non-existent choice
    assert(false)
end

end % of fgnsr_init
