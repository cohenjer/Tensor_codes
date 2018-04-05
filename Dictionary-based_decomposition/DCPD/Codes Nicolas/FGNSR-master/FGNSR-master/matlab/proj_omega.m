function X = proj_omega(C, w, ub, diag_idx, method, verbose)
% X = proj_omega(C, w, ub, diag_idx, method, verbose)
%   -- project matrix C on a certain polyhedron that appears in the 
%      context of an LP formulation for separable NMF.
%
% This function projects the input matrix C \in \R^{n,n} on the polyhedron
%
%   {X in R^{n,n} | w(j) * X(i,j) <= w(i) * X(j,j) forall i,j,
%                   X(j,j) <= ub forall j, X >= 0}
%
% such that the Euclidean distance of each pair of columns
% norm(X(:,j) - C(:,j), 2) is minimum.
%
% Inputs:
%   C       --  Input matrix to project.
%   w       --  Nonnegative vector of weights (def.: all-ones)
%   ub      --  Nonnegative upper bound on the entries X(i,i) (def.: +inf)
%   diag_idx -- Permutation vector for the diagonal entries: diag_idx(j)
%               is the index of the diagonal element of column j
%               default: 1:n
%               (Not available for method 'qp')
%   method  --  string paramter, possible values:
%                 'qp':   solve quadratic programming problem using CPLEX
%                 'thres': slow n*log(n) thresholdung algorithm
%                 'parproj': fast n*log(n) thresholding  algorithm (default)
%                 'clip': Clip column at diagonal element
%                         CAUTION not optimal
%                 'box': Clip all entries to [0,ub], irrespective of w
%   verbose --  true/false switch that controls additional progress output
%
% Output:
%   X       --  the projected matrix
%
% See also: FGNSR


[m,n] = size(C);

if (m < n)
    error('NMFPACK:input', 'C must not be a short fat matrix');
end

if nargin < 2 || isempty(w)
    w = ones(m,1);
end

if nargin < 3 || isempty(ub)
    ub = +inf;
end

if nargin < 4 || isempty(diag_idx)
    diag_idx = 1:n;
end

if nargin < 5 || isempty(method)
    method = 'parproj';
end

if nargin < 6 || isempty(verbose)
    verbose = false;
end

if any(isinf(w)) || any(w < 0.0)
    error('FGNSR:input', 'weights must be finite and non-negative');
end

if (ub < 0)
    error('Upper bound on diagonal elements must be nonnegative');
end

if strcmp(method, 'qp')
    X = qp_cplex(C, w, ub, diag_idx, verbose);
    return;
    % qp_cplex returns whole matrix, unlike other methods. Can exit now.
elseif strcmp(method, 'parproj')
    X = parproj(C, w, ub, diag_idx, verbose);
    return;
elseif strcmp(method, 'box')
    X = clip_to_box(C, ub);
    return;
elseif strcmp(method, 'thres')
    project_col = @project_by_thresholding;
elseif strcmp(method, 'clip')
    project_col = @just_clip;
else
    error('Unknown value for method paramter: %s', method');
end

X = zeros(size(C));
min_weight = 0.0;
sumC = sum(abs(C)); 


for j=1:n
    % Permute diagonal entry to first entry, which is the assumption
    % for the single column projection functions
    z = C(:,j);
    this_diag_idx = diag_idx(j);
    z([1 this_diag_idx]) = z([this_diag_idx 1]);
    w([1 this_diag_idx]) = w([this_diag_idx 1]);
    
    if w(1) > min_weight && sumC(j) > 0
        X(:,j) = project_col(z, w, ub, verbose);
    else
        % The weight for this column is to small, hence no restrictions
        % other than being in the positive orthant and ub.
        z(z<0) = 0.0;
        z(1) = min(z(1), ub);
        X(:,j) = z;
    end
    % Permute back
    X([1 this_diag_idx], j) = X([this_diag_idx 1], j);
    w([1 this_diag_idx]) = w([this_diag_idx 1]);
end
    
end % of main function



function C = qp_cplex(C, weights, ub, diag_idx, verbose)
% Solve QP formuation for the problem using CPLEX
[m,n] = size(C);

E = speye(m,m);

cpx = Cplex();
cpx.Model.Q = 2*speye(m,m);
cpx.Model.lb = zeros(m,1);
cpx.Model.ub = inf * ones(m,1);
cpx.Model.lhs = -inf * ones(m,1);
cpx.Model.rhs = zeros(m,1);
cpx.Model.sense = 'minimize';

cpx.Param.qpmethod.Cur = 2;
cpx.Param.advance.Cur = 0;


if ~verbose
    cpx.DisplayFunc = [];
end

for j=1:n
    pj = diag_idx(j);
    c = C(:,j);
    Aoffset = sparse(1:m,pj*ones(m,1),weights,m,m);
    cpx.Model.obj = -2*c;
    cpx.Model.A =  weights(pj) * E - Aoffset;
    cpx.Model.ub(pj) = ub;
    
    cpx.solve();
    if verbose
        objval = cpx.Solution.objval;
        fprintf('Solution time: %.1ds\n', cpx.Solution.time);
        fprintf('Optimal objval: % .4f\nSquared l2dist: % .4f\n', ...
            objval, objval + sum(c.^2));
    end
    C(:,j) = cpx.Solution.x;
    
    % Clean up objective
    cpx.Model.obj = [];
    cpx.Model.ub(pj) = +inf;
end

end % of qp_cplex


function x = project_by_thresholding(z, w, ub, verbose)

assert(w(1) > 0);

% Remove entries corresponding to zero weight.  For these k we have
%   w(1) * z(k) <= 0 * z(1)
% which implies z(k) = 0.  This will be enforced in 'clip_at_x1' after
% the optimal value of x1 has been determined.  For this determinatino of
% x1, such entries corresponding to zero weight have no effect.
mask = (w>0);
z_org = z;
w_org = w;
z = z(mask);
w = w(mask);

n = length(z);

% These are the upper bounds on z w.r.t. w and x(1).
assert(all(w>0));
values = z ./ w;

% Sort components by (weighted) violation of the condition
% w(1)*x(j) <= w(j)*x(1)
[~, p] = sort(values(2:end), 'descend');

% Initial value for the weighted mean.  If z is feasible, this will be the
% final value.
nom = w(1) * z(1);
den = w(1) * w(1);
x1 = max(0, w(1) * nom/den);
x1 = min(x1, ub);

for k=1:n-1
    which = p(k) + 1;
    if w(1) * z(which) <= w(which) * x1
        % Then the next-largest originally violated constraint is satisfied
        % by the current choice of x1, exit.
        break;
    end
    
    % Find a better x1 such that the cost of satisfying all violated
    % constraints is minimum.
    nom = nom + w(which) * z(which);
    den = den + w(which)^2;
    x1 = max(0,w(1) * nom/den);
    x1 = min(x1, ub);
end

% We determined the optimal value of x(1); given that value, the optimal
% Euclidean projection is found simply by clipping components outside the
% rectangle defined by x(1) to the boundary.
x = clip_at_x1(z_org, w_org, x1);

if verbose
    fprintf('Proj. has %d tight upper bounds implied by x1=%.3f\n',k-1,x1);
end


end % of function proj_threshold

function x = just_clip(z, w, ub, verbose) %#ok<INUSD>
% Project column in a non-optimal fashion:  just clip at x1

x1 = max(0, min(z(1), ub));

x = clip_at_x1(z, w, x1);

end

function x = clip_at_x1(z, w, x1)
% Given point z and a value for x(1), return optimal projected
% point x satisfying box constraints implied by weights w and x(1).

n = length(z);
x = zeros(size(z));
x(1) = x1;

for k=2:n
    x(k) = min(max(0,z(k)), x(1) * w(k) / w(1));
end

end

function X= clip_to_box(C, ub)
% Clip all entries to the interval [0, ub]
assert(isfinite(ub));
X = C;
X(X<0) = 0.0;
X(X>ub) = ub;
end
