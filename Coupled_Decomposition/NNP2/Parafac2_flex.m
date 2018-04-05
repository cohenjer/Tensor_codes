function  [A_aux,C_aux,D_aux,P,C_s,Y_est,err_als_vec] = ...
          Parafac2_flex(Y,iter,A_0,C_0,D_0,opts,P_0)
%-------------------------------------------------------------------------%
% [A_aux,C_aux,D_aux,P,C_s,Y_est,err_als_vec] = ...
%           Parafac2_flex(Y,iter,A_0,C_0,D_0)
%
% Proprocessed ALS algorithm for Coupled CP decomposition for N tensors. Factors Ck 
% are coupled through estimated factors Pk and with coupling error. Factors Ak are coupled
% directly without coupling error.
%
%      cost = sum_k( || T(:,:,k) - A*D_k*C_k' ||^2_F + mu_k ||C_k' - C_s'*P_k' ||^2_F )
%               where D_k is a diagonal matrix stacked vertically in D_aux.
%
% Inputs:
% 
% - Y           : cell of data tensors;
% - iter        : maximum number of iterates;
% - A_0,...     : cell of initial factors;
% - opts        : opts.cstr={'NN','NN','NN'} for nonnegativity constraints;
%                 opts.increase = multiplicative increase of mu;
%                 opts.mu_lvl = amount of regularization;
% - P_0         : cell of initital transformation matrices for the factors
%                 C;
% 
%                 
%
% Outputs:
% 
% - A,C,D       : cell of estimated factors;
% - C_s         : coupling core matrix Ck = Pk*C_s;
% - P           : cell of estimated coupling matrices;
% - map_als_vec : objective function.
%
% List of updates                 -     03/10/2017   J.E.Cohen
%                                       Creation of the file
%                                 -     05/04/2018   J.E.Cohen
%                                       Normalisation, cleanning up
% TODO: - Allow various noise levels for each slices (not mandatory)
%       - Default initialization depending on constraint set.
%       - Slices with different sizes. (easy)
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%

fprintf('\n Flexible Parafac2 for easy constraints \n')

%-------------------------------Parameters--------------------------------%
dim       =     size(Y);
N         =     dim(3);

if nargin<3
    % Ask for the rank
    prompt ='What is the rank of the decomposition ? (type desired integer) ';
    R = input(prompt);
    
    % Random initial factor (NN by default)    
    A_0 = rand(dim(1),R);
    C_0 = rand(dim(2),R);
    D_0 = rand(dim(3),R);
    
    % Options
    opts.cstr     = {'NN','NN','NN'};
    opts.increase = 1.05;
    opts.mu_lvl   = 0.2;
    opts.max_mu   = 10^12;
    for i=1:N
        P_0{i}       =     eye(dim(2),R); % initialize with SVD
    end
elseif nargin<6
    % Tensor estimated rank
    R         =     size(A_0,2);
    
    % Default options
    opts.cstr     = {'NN','NN','NN'};
    opts.increase = 1.05;
    opts.mu_lvl   = 0.2;
    opts.max_mu   = 10^12;
    for i=1:N
        P_0{i}       =     eye(dim(2),R); % initialize with SVD
    end
elseif nargin<7
    % Tensor estimated rank
    R         =     size(A_0,2);
    
    for i=1:N
        P_0{i}       =     eye(dim(2),R);
    end
end

%err_tol = 10^-4; % Error minimal value for stopping criterion
dif_tol = 10^-4; % Relative error minimal value for stopping algorithm
iter_err = iter/10; % number of iterations before computing error
itermax_NN = 100; % Maximal number of nnls sub-iterations

% Storing initial condition
A_aux     =     A_0;

if iscell(C_0)
    C_aux     =     C_0;
else
    C_aux     =     cell(1,N);
    for i=1:N % Only if C_0 input is a matrix.
        C_aux{i}  =     C_0; % Coupled mode, in the preprocessed space
    end
end
D_aux     =     D_0; % Stacked mode, scaling ratios
P         =     P_0;

err_als_vec(1) = 10000; err_als_vec(2) = 1000;

% Slices definition
M     =     cell(1,N);
for i=1:N
M{i}  =     Y(:,:,i);
end

Ynorm     =     sum(Y(:).^2);


% Initial error
mu = ones(1,N);
err_temp = 0;
for i=1:N
    calc_temp = norm(M{i}-A_aux*diag(D_aux(i,:))*C_aux{i}','fro')^2;
    err_temp  = err_temp + calc_temp;
    % setting initial mu to low but not negligeable value
    mu(i)     =     calc_temp/norm(C_aux{i},'fro')^2*(10^-1); 
end
err_init     =     err_temp/Ynorm;
fprintf('-------Initial fit Error %d ----------\n',err_init)
fprintf('----- Mean Initial mu value %d -------\n\n',mean(mu))
    

%--------------------- Outer Loop for global iterations ------------------%
k=0;
while k<iter && (abs(err_als_vec(end)-err_als_vec(end-1))/err_als_vec(end))>dif_tol
%-----------------------------------------------------%
k=k+1; % iteration increment

% if last iteration
if k==iter
    fprintf('Last Iteration\n')
    itermax_NN = 500;
end
    
%-------------------------------------------------------------------------%
% Coupling increase
mu = min(mu*opts.increase,opts.max_mu);
%------------------------- Coupling Estimation ---------------------------%

% C_s estimate (1)
% Av_C = 0;
% for i=1:N
%     Av_C = mu(i)*C_aux{i}'*C_aux{i}+Av_C;
% end
% C_s = sqrtm(Av_C/sum(mu));

for l=1:5 % This is fast, so let's do it more times
% C_s estimate (2)
Av_C = 0;
for i=1:N
    Av_C = Av_C + mu(i)*P{i}'*C_aux{i};
end
C_s = Av_C/sum(mu);
C_s = col_norm(C_s);
    
% Pk estimates
for i=1:N
    [u,~,v] = svd(C_aux{i}*C_s','econ');
    P{i}    = u*v';
end
end

%--------------------------ALS Steps--------------------------------------%
    
%-------------------------- NALS one iteration-----------------------------%
       
    % A update
    temp = 0;
    temp_inv = 0;
    for i=1:N
        temp = temp + M{i}*C_aux{i}*diag(D_aux(i,:));
        temp_inv = temp_inv + diag(D_aux(i,:))*C_aux{i}'*C_aux{i}*diag(D_aux(i,:));
    end
    if strcmp(opts.cstr{1},'NN')
    A_aux = nnlsHALSupdt(temp',temp_inv',A_aux',itermax_NN)';
    else
    A_aux = temp/temp_inv;
    end
    [A_aux,D_aux] = col_norm(A_aux,D_aux);
    
    
    % C update (funny mode) 
    for i=1:N
    if strcmp(opts.cstr{2},'NN')
        C_aux{i} = nnlsHALSupdtv2(M{i},...
            A_aux*diag(D_aux(i,:)),C_aux{i}',itermax_NN,C_s'*P{i}',mu(i))';
    else
        C_aux{i} = (M{i}'*A_aux*diag(D_aux(i,:))+mu(i)*P{i}*C_s)/(...
         diag(D_aux(i,:))*(A_aux'*A_aux)*diag(D_aux(i,:)) + mu(i)*eye(R));
    end
    end
    
    % D factors update (Slice by Slice) 
    for i=1:N    
    if strcmp(opts.cstr{3},'NN')
        D_aux(i,:) = nnlsHALSupdt(M{i}(:),kr(C_aux{i},A_aux),D_aux(i,:)',itermax_NN)';
    else
        D_aux(i,:) = diag((A_aux\M{i})/(C_aux{i}'))';
        %D_aux(i,:) = (pinv((B{k}'*B{k}).*(A'*A))*(kr(B{k},A)'*M{k}(:)))'; 
    end
    end
   % normalization ??
   %D_aux = col_norm(D_aux);
    

    if mod(k,iter_err)==1 
        
            % Error computation (to be placed later in the code after debugging)
    err_temp = 0;
    err_fit  = 0;
    temp_i_r = zeros(1,N);
    for i=1:N
        fit_i =  norm(M{i}-A_aux*diag(D_aux(i,:))*C_aux{i}','fro')^2;
        temp_i = mu(i)*norm(C_aux{i} - P{i}*C_s,'fro')^2;
        temp_i_r(i) = temp_i/norm(C_aux{i},'fro')^2/mu(i);
        % second iteration to have an idea of the amplitudes of the terms in the cost function
        if k==1 
            mu(i) = opts.mu_lvl*fit_i/temp_i;
            %fprintf('new value of mu: %d\n',mu(i))
        end
        err_fit  = err_fit + fit_i;
        err_temp = err_temp + temp_i;
    end
    err_als_vec(floor((k-1)/iter_err)+1)    =     (err_fit+err_temp)/Ynorm;
        
    fprintf('Iteration %g\n Cost function          %d\n',k,err_als_vec(floor((k-1)/iter_err)+1))
    fprintf(' Relative fitting error %d\n',err_fit/Ynorm)
    fprintf(' Avg. Rel. coupl. error %d\n',mean(temp_i_r))
    fprintf(' Max. Rel. coupl. error %d\n',max(temp_i_r))
    fprintf(' Mean mu value          %d\n',mean(mu))
    end
    
    % Estimated tensor
    
    Y_est     =     zeros(size(Y));
    
    for i=1:N
        Y_est(:,:,i)     =       A_aux*diag(D_aux(i,:))*C_aux{i}';
    end
  
end % end while

end % end PARAFAC2_flex

function V = nnlsHALSupdt(M,U,V,maxiter) 
    % Credits to Nicolas Gillis    
% See N. Gillis and F. Glineur, Accelerated Multiplicative Updates and 
% Hierarchical ALS Algorithms for Nonnegative Matrix Factorization, 
% Neural Computation 24 (4): 1085-1105, 2012.

if nargin <= 3
    maxiter = 500;
end
[m,n] = size(M); 
[m,r] = size(U); 
UtU = U'*U; 
UtM = U'*M;

if nargin <= 2 || isempty(V) 
    V = U\M; % Least Squares
    V = max(V,0); 
    alpha = sum(sum( (U'*M).*V ))./sum(sum( (U'*U).*(V*V'))); 
    V = alpha*V; 
end

delta = 1e-6; % Stopping condition depending on evolution of the iterate V: 
              % Stop if ||V^{k}-V^{k+1}||_F <= delta * ||V^{0}-V^{1}||_F 
              % where V^{k} is the kth iterate. 
eps0 = 0; cnt = 1; eps = 1; 
while eps >= (delta)^2*eps0 && cnt <= maxiter %Maximum number of iterations
    nodelta = 0; if cnt == 1, eit3 = cputime; end
        for k = 1 : r
            deltaV = max((UtM(k,:)-UtU(k,:)*V)/UtU(k,k),-V(k,:));
            V(k,:) = V(k,:) + deltaV;
            nodelta = nodelta + deltaV*deltaV'; % used to compute norm(V0-V,'fro')^2;
            if V(k,:) == 0, V(k,:) = 1e-16*max(V(:)); end % safety procedure
        end
    if cnt == 1
        eps0 = nodelta; 
    end
    eps = nodelta; 
    cnt = cnt + 1; 
end

end % of function nnlsHALSupdt
  
function V = nnlsHALSupdtv2(M,U,V,maxiter,Vtarget,lambda)
    % Credits to Nicolas Gillis
if nargin <= 3
    maxiter = 500;
end
[~,n] = size(M); 
[~,r] = size(U); 
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
            if V(k,:) == 0 
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
end % of function nnlsHALSupdtv2

function [ M,N ] = col_norm( M,N )
[m,n] = size(M);
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
end % end of col_norm

function X = kr(U,varargin)
%   Version: 21/10/10
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)

if ~iscell(U), U = [U varargin]; end
K = size(U{1},2);
if any(cellfun('size',U,2)-K)
    error('kr:ColumnMismatch', ...
          'Input matrices must have the same number of columns.');
end
J = size(U{end},1);
X = reshape(U{end},[J 1 K]);
for n = length(U)-1:-1:1
    I = size(U{n},1);
    A = reshape(U{n},[1 I K]);
    X = reshape(bsxfun(@times,A,X),[I*J 1 K]);
    J = I*J;
end
X = reshape(X,[size(X,1) K]);
end % end kr