function [X, K, fval, mmx_err, fval_hist, diag_hist] = fgnsr(M, r,varargin)
% [X, K, fval, mmx_err, fval_hist, diag_hist] = fgnsr(M, r, p, ...)
%
% Fast gradient method for nonnegative sparse regression with self dictionary.
% For a given input matrix M, it solves the optimization problem
%
%   min_{X in Omega} ||M-MX||_F^2 + mu * p^T diag(X),
%
% where
%
% Omega = {X in R^{n,n} | w(j) * X(i,j) <= w(i) * X(j,j) forall i,j,
%                                X(j,j) <= ub forall j, X >= 0}
%
% for some weights w(i), which are typically the l1-norms of the columns of
% the input matrix M.
%
% Details can be found at
%
%     https://arxiv.org/TODO
%
% Inputs:
%   M  -- the input matrix
%   r  -- number of separable columns to extract
%   p  -- objective coefficients for the trace of X
%         (optional, default: ones(n,1) + .01 * (rand(n,1) - .5) )
%
% Further arguments can be passed in ('identifier', value) pairs:
%
%   'maxiter' -- integer iteration limit (default: size(M,2))
%
%   'init'    -- initialization method used.  Possible values:
%       0: use random K for mu, set X = 0
%       1: use FastSepNMF K for mu, set X = 0 (default)
%       2: use FastSepNMP K + KKT heuristic for mu, set X=0
%
%   'verbose' -- true/false switch to control output (default: true)
%
%   'log_interval' -- print every n-th iteration log line (default: 25)
%
%   'mu'      -- Value for scaling factor mu, overrides the mu implied by
%     'init' option (default: Use mu from 'init').
%
%   'dynamic_mu' -- integer switch to control the Lagrangian variable mu:
%        -1:  do not adjust mu
%        0:   adjust mu upon restart (default)
%        k>0: adjust mu every k iterations
%
%   'target' -- Define target quantity for dynamic mu adjustments
%        'residual' -- aim for certain M-MX residual
%        'trace'    -- aim for certan tr(X) value (default)
%
%   'trX_target' -- targeted value of Tr(X) for mu adjustments
%        (default: r/2+1)
%
%   'res_target' -- targeted value of |M-MX|_F
%        (default: estimate by SVD / FastSepNMF)
%
%   'restart' -- integer, restart after every k iterations,
%                   0 means no restart (default)
%
%   'subset'  -- integer array of column indices of M to which the solution
%                should be restricted to (default: 1:n)
%
%   'pspread' -- spread for the coefficients in the objective vector p:
%          ones(n,1) + pspread * (rand(n,1) - .5)  (default: 0.01)
%
%   'normalize' -- Normalize the columns of the input matrix by the
%                column sums. (true/false switch, default: false)
%
%   'proj'    -- projection type, possible values
%        'thres'  : Slow thresholding algorithm
%        'qp'     : Very slow explicit QP based projection (requires CPLEX)
%        'parproj': Very fast thresholding algorithm (default)
%        'clip'   : Extremely fast but NON OPTIMAL clipping heuristic
%        'box'    : Clip all entries only to the box [0,1]
%                   CAUTION: this makes only sense if the colunms of M
%                   are normalized.
%
%   'col_weights' -- column weights for the projection
%                    default: sum(abs(M))
%
% The iteration log reads as follows:
%
% iter  chg_diag        fval  |M-MX|_F            mu   tr(X)  nnz   sep
%    0  0.00e+00  2.4145e+00  4.8290e+00  3.6874e-04    0.00    0  0.00
%   50  1.06e-02  9.7413e-02  6.0478e-02  2.2534e-04    2.98   51  0.03
%  100  5.81e-03  1.0637e-01  5.1587e-02  2.6822e-04    3.00   36  0.05
%  150  4.02e-03  9.6126e-02  5.0995e-02  2.3844e-04    2.96   29  0.06
%  199  2.52e-03  9.5896e-02  4.6569e-02  2.3962e-04    3.03   21  0.08
%
% Meaning of the fields:
%
%   iter     -- iteration count
%   chg_diag -- 1-norm of changes to the diagonal of X from one iteration
%               to the next, relative to current 1-norm of diag(X)
%   fval     -- current value of  0.5 * ||M-MX||_F^2 + p^T diag(X)
%   |M-MX|_F^2 -- absolute error
%   trX      -- trace(X)
%   nnz      -- nnz(diag(X))
%   sep      -- relative gap from the r-th largest diagonal value of X
%               to the (r+1)st largest value

[m,n] = size(M); 

ip = inputParser;

ip.addOptional('p', [], @(x) length(x)==n);

ip.addParamValue('maxiter', n, @(x) 0<=x);
ip.addParamValue('init', 1, @(x) 0<=x && x<=5);
ip.addParamValue('verbose', true, @(x) islogical(x));
ip.addParamValue('log_interval', 25, @(x) 1 <= x);
ip.addParamValue('gradient', 'full', ...
    @(x) strcmp(x, 'full') || strcmp(x, 'sparse'));
ip.addParamValue('mu', [], @(x) isscalar(x) && x>=0);
ip.addParamValue('dynamic_mu', 0, @(x) x >= -1);
ip.addParamValue('target', 'trace', ...
    @(x) strcmp(x,'trace') || strcmp(x, 'residual'));
ip.addParamValue('trX_target', r/2+1, @(x) x >= 0);
ip.addParamValue('res_target', [], @(x) isempty(x) || x>= 0);
ip.addParamValue('restart', 0, @(x) x >= 0);
ip.addParamValue('subset', []);
ip.addParamValue('normalize', false, @(x) islogical(x));
ip.addParamValue('pspread', 0.01, @(x) x >= 0 && x <= 1.0);
ip.addParamValue('proj', 'parproj', @(x) strcmp(x, 'thres') || ...
    strcmp(x, 'qp') || strcmp(x, 'clip') || strcmp(x, 'parproj') || ...
    strcmp(x, 'box'));
ip.addParamValue('col_weights', sum(abs(M)), @(x) all(x >= 0) && ~isempty(x));

ip.parse(varargin{:});
verbose = ip.Results.verbose;

p = ip.Results.p;
if isempty(p)
    pspread = ip.Results.pspread;
    if verbose
        fprintf ('Using coefficient spread of %f for p\n', pspread);
    end
    p = ones(n,1) + pspread * (rand(n,1) - .5);
end

opt = ip.Results;

maxiter = ip.Results.maxiter;
init = ip.Results.init;
log_interval = ip.Results.log_interval;
dynamic_mu = ip.Results.dynamic_mu;
trX_target = ip.Results.trX_target;
restart = ip.Results.restart;
subset = ip.Results.subset;
proj_type = ip.Results.proj;
gradient = ip.Results.gradient;


if opt.normalize
    inv_col_sums = 1.0./(sum(abs(M)) + 2*eps);
    M = M * spdiags(inv_col_sums', 0, size(M,2), size(M,2));
end

if isempty(subset)
    Ms = M;
else
    % Put subset columns to leading part of M
    subset = subset(:)';
    mask = true(1,n);
    mask(subset) = false;
    tail = find(mask==true);
    head = subset;

    colperm = [head, tail];
    M = M(:, colperm);
    Ms = M(:, 1:length(subset));
    p = p(subset);
end

[ms,ns] = size(Ms);

rold = r;
J = FastSepNMF(Ms, r, 0);
if length(J) < r
    warning('NN rank is at most %d\n', length(J)); %#ok<*WNTAG>
    r = length(J);
end

if trX_target > r
    warning('trX_target larger than (poss. reduced) r');
    trX_target = r;
end


if nargout >= 4
    diag_hist = zeros(maxiter+1,ns);
else
    diag_hist = [];
end

if r == 0
    % The zero matrix was input...
    X = zeros(ms ,ns);
    K = [];
    fval = 0;
    mmx_err = 0;
    fval_hist = 0;
    return;
end

iter_log = iteration_log_factory(verbose, @print_header, @print_line, ...
    log_interval);

init_X = 0;
init_K = init;

[X, ~, mu] = fgnsr_init(Ms, r, p, init_K, init_X);
% Extend X to other columns of M
X(:, ns+1 : n) = 0.0;

if ~isempty(ip.Results.mu)
    mu = ip.Results.mu;
end


target = ip.Results.target;
if strcmp(target, 'residual')
    if isempty(ip.Results.res_target)
        svdvals = svds(Ms, r);
        svd_est = sqrt(norm(Ms, 'fro')^2 - sum(svdvals.^2));
        res_target = 0.5 * ( svd_est + sqrt(spa_err_squared) );
    else
        res_target = ip.Results.res_target;
    end
end

if strcmp(gradient, 'sparse')
    updtset = false(ns, 1);
    max_updt_size = 5*r;
    candsize = min(ns, 2*r);
    updt_interval = min(25, floor(maxiter / (2 * r)));
    % While we do not have enough columns in the update set, we set mu
    % to zero and disallow any mu adjustments.  As soon there are at least
    % r columns in the update set, we switch to the regular mu and update
    % it if the the other settings say so.
    mu_phase = false;
    mu_value = mu;
    mu = 0;
else
    updtset = true(ns,1);
    mu_phase = true;
end

% l-norm of the columns of M 
w = ip.Results.col_weights;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See Constant step scheme, II - in "Introductory Lectures on Convex
% Programming, Volume I: Basic course", Yurii Nesterov, 2004 -p.90 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precompute some parts of the gradient
MstM =  Ms' * M; 
MstMs = Ms' * Ms;

pmat = zeros(ns,n);
for k=1:ns
    pmat(k,k) = p(k);
end

% Estimate Lipschitz constant L 
L = svds(Ms,1)^2;
step_len = 1/L;

% Preallocate some array
fval_hist = -1 * ones(maxiter+1,1);
if ~isempty(diag_hist)
    diag_hist(1, :) = diag(X);
end

% Initialization of iteration paramters
Y = X;
alpha1 = 0.05;
num_iter = 0;
was_restarted = false;

[sep, mmx_err, ptx, fval, diag_chg, trX, nnzX, upd] = ...
    iteration_stats(X, zeros(size(X)), M, r, p, mu, updtset);

fval_hist(1) = fval;
if ~isempty(diag_hist)
    diag_hist(1, :) = diag(X);
end

if verbose == true
    fprintf('M has %d rows, %d columns %d nonzeros (density %.3f)\n', ...
        m, n, nnz(M), nnz(M)/(m*n));
    fprintf('Ms has %d rows, %d columns %d nonzeros (density %.3f)\n', ...
        ms, ns, nnz(Ms), nnz(Ms)/(ms*ns));

    fprintf('Ms''M has %d NZ, density %.3f, fill %.2f\n', ...
        nnz(MstM), nnz(MstM)/(ns*n), nnz(MstM)/nnz(Ms));
    fprintf('Ms''Ms has %d NZ, density %.3f, fill %.2f\n', ...
        nnz(MstMs), nnz(MstMs)/(ns*ns), nnz(MstMs)/nnz(Ms));
    
    fprintf('Estimated Lipschitz constant is %.2e\n', L);

    if strcmp(target, 'trace')
        fprintf('Target trace value is %.2f\n', trX_target);
    else
        fprintf('Target residual value is %.3e\n', res_target);
    end
    
    fprintf('Projection method is %s\n', proj_type);
end


iter_log(num_iter, diag_chg, fval, mmx_err, ptx, mu, trX, nnzX, sep, ...
    was_restarted, upd);

% Needed for dynamic mu adjustments
sigma_dir = 0;
sigma = 0.1;

force_restart = false;
Xold = X;

if strcmp(gradient, 'sparse')
    new_row_idx = choose_one_row(updtset, candsize, X, MstM, ...
                MstMs, mu, pmat, w, proj_type);
    updtset(new_row_idx) = true;
end
while num_iter < maxiter
    num_iter = num_iter + 1;

    % alpha2 is solution of alpha2^2 = (1-alpha2)*alpha1^2
    alpha2 = ( sqrt(alpha1^4 + 4*alpha1^2 ) - alpha1^2) / (2); 
    beta = alpha1*(1 - alpha1)/(alpha1^2 + alpha2); 

    if strcmp(gradient, 'full')
        Xold = X;

        gradY = -MstM + MstMs * Y + mu * pmat;

        % Update iterates
        X = proj_omega((Y - step_len * gradY)', w', 1, [], proj_type)'; 
        Y = X + beta * (X - Xold); 


    else
        assert(strcmp(gradient, 'sparse'));
        % Update only the rows designated by updtset
        Xold(updtset,:) = X(updtset,:);
        gradY = - MstM(updtset,:) + MstMs(updtset,updtset) * Y(updtset,:) + mu * pmat(updtset,:);
        X(updtset,:) = proj_omega( (Y(updtset,:) - step_len * gradY)', w', 1, find(updtset), proj_type)';
        Y(updtset,:) = X(updtset, :) + beta * (X(updtset,:) - Xold(updtset,:));

        if mod(num_iter, updt_interval) == 0
            x = diag(X);
            delcand = find((x > 0) & (x < 1e-4) & (updtset==true));
            if ~isempty(delcand)
                X(delcand,:) = 0.0;
                updtset(delcand) = false;
                force_restart = true;
            end
        end
        
        if mod(num_iter, updt_interval) == 0 && upd < max_updt_size
            % Increase the set set of rows that are updated
            
            new_row_idx = choose_one_row(updtset, candsize, X, MstM, ...
                MstMs, mu, pmat, w, proj_type);
            
            % We may not find any useful new row
            if ~isempty(new_row_idx)
                updtset(new_row_idx) = true;
                % If there is a new row, forget history
                force_restart = true;
            end
            
            if nnz(updtset) >= r && mu_phase == false
                mu_phase = true;
                mu = mu_value;
            end
        end
    end

    % Iteration upkeep
    fval_old = fval;
    [sep, mmx_err, ptx, fval, diag_chg, trX, nnzX, upd] = ...
        iteration_stats(X, Xold, M, r, p, mu, updtset);

    if restart == 0 && ~was_restarted && fval_old < fval
        force_restart = true;
    end
    
    if (restart > 0 && mod(num_iter, restart) == 0) || force_restart
        force_restart = false;
        Y = X;
        alpha1 = 0.05;
        was_restarted = true;
    else
        alpha1 = alpha2;
        was_restarted = false;
    end
    
    
    % Save iteration history
    fval_hist(num_iter + 1) = fval;
    if ~isempty(diag_hist)
        diag_hist(num_iter + 1, :) = diag(X);
    end
    
    iter_log(num_iter, diag_chg, fval, mmx_err, ptx, mu, trX, nnzX, sep,...
        was_restarted, upd);

    if mu_phase && ( (dynamic_mu == 0 && was_restarted) || ...
            (dynamic_mu > 0 && mod(num_iter, dynamic_mu) == 0) )
        
        % Update logic for chasing tr(X) or mmx_err values.
        % Note both are updated by the same update logic, only the first
        % two arguments are swapped to reflect the reciprocal effect
        % between the two quantities in the mu adjustment.
        if strcmp(target, 'residual')
            [mu, sigma, sigma_dir] = update_mu(mmx_err, res_target, mu, sigma, ...
                sigma_dir);
            
        else
            [mu, sigma, sigma_dir] = update_mu(trX_target, trX, mu, sigma, ...
                sigma_dir);
        end
    end

    if diag_chg < 1e-8 && num_iter > 1
%        break;
    end
    
end

% Force a log line for the last feeded iteration data ("flush")
iter_log();

[~, K] = sort(diag(X), 'descend');
K = K(1:rold);
if ~isempty(subset)
    K = subset(K);
end

K = reshape(K, 1, length(K));

fval_hist = fval_hist(1:num_iter+1);
if ~isempty(diag_hist)
    diag_hist = diag_hist(1:num_iter+1, :);
end

end % of fgnsr

function print_header
fprintf('%s%4s  %8s  %9s  %9s  %9s  %9s  %6s  %3s  %4s  %3s\n', ...
    ' ', 'iter', 'chg_diag', 'fval', '|M-MX|_F', 'p^T diagX', 'mu', ...
    'tr(X)', 'nnz', 'sep', 'upd');
end

function print_line(iter, chg_diag, fval, mmx_err, ptx, mu, trX, nnzX, ...
    sep, restart, upd)
if restart == true
    restart = '+';
else
    restart = ' ';
end

fprintf('%s%4d  %.2e  %.3e  %.3e  %.3e  %.3e  %6.2f  %3d  %4.2f  %3d\n', ...
    restart, iter, chg_diag, fval, mmx_err, ptx, mu, trX, nnzX, sep, upd);
end

function fhandle = iteration_log_factory(verbose, header_func, ...
    iter_func, log_interval)
if nargin < 4 || isempty(log_interval)
    log_interval = 10;
end

is_flushed = false;
cache = {};

header_interval = 24;
lines_printed = 0;
num_calls = 0;

    function display_function(varargin)
        % If no parameters are passed, it means that the cache is to be
        % flused
        if isempty(varargin)
            flush = true;
        else
            flush = false;
        end

        if flush == true
            data = cache;
        else
            data = varargin;
        end
        
        if num_calls == log_interval
            num_calls = 0;
        end
        
        if num_calls == 0 || (flush == true && ~is_flushed)
            if lines_printed == header_interval
                lines_printed = 0;
            end
            if lines_printed == 0
                header_func();
            end
            
            iter_func(data{:})
            lines_printed = lines_printed + 1;
            is_flushed = true;
        else
            is_flushed = false;
        end
        
        cache = varargin;
        num_calls = num_calls + 1;
    end

    function empty_function(varargin)
    end

if verbose == true
    fhandle = @display_function;
else
    fhandle = @empty_function;
end
end

function [sep, mmx_err, ptx, fval, diag_chg, trX, nnzX, upd] = ...
        iteration_stats(X, Xold, M, r, p, mu, updtset)
    [m,~] = size(X);
    values = sort(diag(X), 'descend');

    if isempty(values)
        sep = 0;
    elseif values(r) == 0
        sep = 0;
    elseif r==m
        sep = 1;
    else
        sep = (values(r) - values(r+1)) / values(r);
    end
    
    mmx_err = norm(M - M(:,updtset) * X(updtset,:), 'fro');
    
    if m > 1
        ptx = p'*diag(X);
    else
        ptx = p * X(1,1);
    end
    muptx = mu * ptx;
    fval = 0.5 * mmx_err^2 + muptx;
    normX = norm(diag(X), 1);
    if normX == 0.0
        % This should only happen at initialization
        diag_chg = norm(diag(Xold), 1);
    else
        diag_chg = norm(diag(X)-diag(Xold), 1) / norm(diag(X),1);
    end
    trX = sum(diag(X));
    nnzX = nnz(diag(X));
    upd = nnz(updtset);
end

function [mu, sigma, sigma_dir] = update_mu(target, trX, mu, sigma, sigma_dir)
    if trX < target
        if sigma_dir < 0
            sigma = sigma / 2.0;
        end
        mu = (1.0 - sigma) * mu;
        sigma_dir = 1;
    else
        if sigma_dir > 0
            sigma = sigma / 2.0;
        end
        mu = (1.0 + sigma) * mu;
        sigma_dir = -1;
    end

end

function new_row_idx = choose_one_row(updtset, candsize, X, MstM, MstMs, mu, pmat, w, ptype)

J_out = find(updtset == false);
J_in = find(updtset == true);

MJ = MstMs(:, J_in);
XJ = X(J_in, :);

% Compute the diagonal of the gradient for inactive rows
grad_diag = - diag(MstM) + mu * diag(pmat);
grad_diag = grad_diag(J_out);
for k=1:length(grad_diag);
    grad_diag(k) = grad_diag(k) + MJ(k,:) * XJ(:,k);%MstMs(k,:) * X(:,k);
end

% Select some candidates of large *negative* diagonal value
[~, Is] = sort(grad_diag);
candsize = min(candsize, length(Is));
cand_idx = Is(1:candsize);
cand_idx = J_out(cand_idx);

% Evaluate the gradient w.r.t. candidate rows
cand_gradX = -MstM(cand_idx,:) + MstMs(cand_idx, updtset) * X(updtset,:) + ...
    mu * pmat(cand_idx,:);

% Project the gradient (the rows of X(cand_idx,:) are zero)
% Also, no step length multiplier necessary as the projection
% is scale invariant.
cand_X = proj_omega( -cand_gradX', w', 1, cand_idx, ptype)';
diag_vals = diag(cand_X(:, cand_idx));

% Only rows *not* projected on zeros are useful
if nnz(diag_vals > 0) == 0
    new_row_idx = [];
    return;
end

% Else: pick the largest one
row_sums = sum(cand_X, 2);
rating = diag_vals .* (row_sums - diag_vals);
[~, which] = max(rating);
new_row_idx = cand_idx(which);

end

