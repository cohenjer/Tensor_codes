% This codes tries to solve NMF with dictionary: 
% given a nonnegative matrix X (p x n), a dictionary D (p x d) and 
% an integer r, find a nonnegative matrix V (r x n) and an index set K 
% with r elements such that D(:,K) V is close to X. 
% 
% min_{K,V>=0} || X - D(:,K) V ||_F^2 such that K in {1,2,...,d}, |K|=r. 
% 
% Our algorithm introduces an auxiliary variable U, close to D(:,K). 
% See Cohen and Gillis, Dictionary-based nonnegative matrix factorization,
% 2017. 
% 
% *** Input ****
% X       :  nonnegative matrix to approximate
% D       : dictionary 
% r       : number of atoms to select in the dictionary
% maxiter : maximum number of iterations (default = 30)
% init  : initialization for the index set K of selected atoms.
%         =1: SPA, 
%         =2: VCA, 
%         -For SNPA and H2NMF, you need to download the codes separately-
%         - available from https://sites.google.com/site/nicolasgillis/ -
%         =3: SNPA, 
%         =4: H2NMF, 
%         ---
%         =5: random initialization, 
%         =vector of r indices: K = vector. 
% 
% *** Output ****
% (K,V)   :  X approx D(:,K)V, with V >= 0, |K| = r 
%  U      :  U approx D(:,K) and minmizes ||X-UV||_F
% re      : evolution of the relative error in percent
% eK      : number of changes in K at each iteration 

function [K,V,U,re,eK] = NMFdico(X,D,r,maxiter,init)

if nargin <= 3 || isempty(maxiter)
    maxiter = 30;
end
if nargin <= 4
    init = 1; % SPA 
end

% Numer of inner iterations in the coordinate descent algorithm 
% to update U and V. 
inneriter = 10; 

% Initialization of index set K 
if length(init) == 1
    if init == 1
        % Initialize with SPA
        K = FastSepNMF(D,r); 
    elseif init == 2
        % Initialize with VCA
        [~, K] = VCA(X,'Endmembers',r); 
    elseif init == 3
        % Initialize with SNPA
        K = SNPA(X,r);  
    elseif init == 4    
        % Initialize with H2NMF
        [~, ~, K] = hierclust2nmf(X,r); 
    elseif init == 5 
        % Random init 
        n = size(D,2); 
        K = randperm(n); 
        K = K(1:r);  
    end
elseif length(init) == r
    K = init;
else 
    error('The parameter ''init'' was not well specified'); 
end
    
% Initialization of U and V 
U = D(:,K); 
V = nnlsHALSupdt(X,U,[],1);  
nX = norm(X,'fro'); 
% Few updares of the NMF algorithm HALSacc 
[U,V] = v2_HALSacc(X,U,V,0.5,0.1,inneriter); 

% Initial value for delta  
delta0 = 0.01 * norm(X-U*V,'fro')^2 / (norm(U-D(:,K),'fro')^2+1e-6); 
delta = delta0; 

% Normalized dico, this will be useful to assign the best atoms to U
p = size(D,1); 
nD = sqrt(sum(D.^2)); 
Dn = D./(repmat(nD,p,1)+1e-6); 

% Main loop 
i = 1; 
eK = 1; 
lastit = 5; 
disp('Display of iteration number and relative error ||X-D(:,K)V||_F/||X||_F*100: ');   
while i <= maxiter 
    Kp = K; 
    Vp = V; 
    
    % Update K 
    for k = 1 : r
        [~,K(k)] = max( Dn'*U(:,k) ); 
    end
    eK(i) = length(setdiff(K,Kp));  
    
    % Update V 
    V = nnlsHALSupdt(X,D(:,K),V,inneriter); 
    % relative error in %  
    re(i) = norm(X-D(:,K)*V,'fro')/nX*100;  
    
    % Update U 
    U = nnlsHALSupdtv2(X',V',U',inneriter,D(:,K)',delta); 
    U = U'; 
    
    % Increase delta 
    if norm(U-D(:,K))/norm(U) > 0.01
        delta = delta*1.5; 
    end
    
    % Display
    fprintf('%2.0f:%2.3f - ',i,re(i));
    if mod(i,5) == 0
        fprintf('\n');
    end
    if i > lastit  &&  sum( eK(i-lastit:i-1) ) == 0  &&  norm(U-D(:,K))/norm(U) < 0.05
        disp('The algorithm has converged.');
        return; 
    end
    i = i+1; 
end
fprintf('\n'); 