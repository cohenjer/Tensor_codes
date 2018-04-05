% Synthetic Data - Experiments on Near-Separable NMF Algorithms
% 
% See Section 4.1. in 
% 
% N. Gillis, Successive Nonnegative Projection Algorithm for Robust 
% Nonnegative Blind Source Separation, arXiv:1310.7529, 2013.  
% 
%
% To reproduce the same experiments as in the paper, use: 
%
% ***************
% Rank-deficient
% ***************
% Dirichlet     : xp = 1, m = 10; r = 20, n = 200, iter = 25, 
%                 delta=logspace(-3,0,100). 
% 
% Middle points : xp = 2, m = 10; r = 20, iter = 25, 
%                 delta=logspace(-3,0,100). 
%
% ***************
% Ill-conditioned
% ***************
%  Dirichlet    : xp = 3, m = r = 20, n = 200, iter = 25, 
%                 delta = logspace(-4,-0.5,100). 
%
% Middle points : xp = 4, m = r = 20, iter = 25, 
%                 delta = logspace(-4, -0.5, 100). 

clear all; clc; close all; 

% Choice of the Experiment: 
xp = 1; 
% 1 => rank-deficient `Dirichlet', 
% 2 => rank-deficient `Middle points', 
% 3 => ill-conditioned `Dirichlet', and 
% 4 => ill-conditioned `Middle points'. 
disp('*************************************************');
if xp == 1
    disp('    Experiment: rank-deficient `Dirichlet'' ');
elseif xp == 2
    disp('    Experiment: rank-deficient `Middle points''');
elseif xp == 3
    disp('    Experiment: ill-conditioned `Dirichlet''');
elseif xp == 4
    disp('    Experiment: ill-conditioned `Middle points''');
end
disp('***********************************************');

% Number of columns of W
r = 10; 
% Number of matrices generated for each noise level
iter = 1;
% Parameter delta (= noise levels)
delta = logspace(-3, -0.5, 10); 
%delta = 1e-3; 

% Dimensions: number of columns and rows of M
if xp == 1 || xp == 2 % Rank-deficient case
    m = 5; % #rows of M, must be smaller than r (for rank deficiency)
elseif xp == 3 || xp == 4 % Ill-conditioned case
    m = 10; % #rows of M, must be larger than r (to have full rank W)
end
% For Middle Points: #col of M = r + nchoosek(r,2), otherwise 
if xp == 1 || xp == 3 % for Dirichlet 
    n = 200; % n is the number of columns of H', with #col of M = n + 2*r. 
end

T = length(delta); 
nalgo = 4; 
results = zeros(nalgo,T);  
resultsree = zeros(nalgo,T);  
timings = zeros(nalgo,1); 
% Construction diagonal matrix for ill-conditioning
param = 3; % kappa(W) = 10^param
alpha = 10^(-param/(r-1)); 
for i = 1 : r
    S(i) = alpha^(i-1);
end

fprintf('Total number of noise levels: %2.0f \n', T)
for i = 1 : T
    for k = 1 : iter
        % Rank-deficient matrix W
        if xp == 1 || xp == 2
            alpha = 0;
            while alpha < 1e-2
                W = rand(m,r); 
                alpha = alphaparam(W); 
            end
        % Ill-conditioned matrix W
        elseif xp == 3 || xp == 4
                W = rand(m,r); 
                [u,s,v] = svds(W,r);
                W = u*diag(S)*v'; 
                W = max(W,0); 
        end
        % Dirichlet experiment
        if xp == 1 || xp == 3
            % H is generated so that the columns of W are repeated twice
            % and the remaining columns are drawn following a Dirichlet 
            % distribution 
            alpha = rand(r,1);  % Parameter of the Dirichlet distribution
            H = [eye(r) eye(r) sample_dirichlet(alpha,n)']; 
            Noise = delta(i)*[randn(m,n+2*r)];  % Noise is Guassian
            M = W*H; 
        % Middle points experiment
        elseif xp == 2 || xp == 4
            % H is generated so that the first r columns of M are the columns 
            % of W, and the remaining ones are on the middle point of any 
            % combination of two columns of W
            H = [eye(r) nchoose2(r)/2]; 
            [r,n] = size(H); 
            M = W*H; 
            Noise = delta(i)*[zeros(m,r) M(:,r+1:end)-repmat(mean(W,2),1,n-r)]; 
        end
        % Input matrix 
        M = W*H + Noise; 
        nM = norm(M,'fro'); 
        % Running the different algorithms
        for algoix = 1 : nalgo
            e = cputime; 
            if algoix == 1 % SPA
                K = FastSepNMF(M,r); 
            elseif algoix == 2 % SNPA
                K = SNPA(M,r); 
            elseif algoix == 3 % XRAY
                K = FastConicalHull(M,r)'; 
            elseif algoix == 4 % dicoNMF
                cd('C:\Users\nicolas\Dropbox\Jérémy-ERC\Papier Dictionnaires\EUSIPCO\matlab'); 
                K = NMFdico(M,M,r,30); 
                cd('C:\Users\nicolas\Dropbox\Jérémy-ERC\Papier Dictionnaires\EUSIPCO\matlab\snpa_v3'); 
            end
            % Compute performance and running time 
            timings(algoix) = timings(algoix) + cputime-e; 
            Ks = K; 
            if xp == 1 || xp == 3
                K1 = K(K <= r); K2 = K(K>=r+1)-r; K = [K1 K2]; 
            end
            rcur = sum(unique(K) <= r); 
            results(algoix,i) = results(algoix,i) + rcur/iter/r; 
            V = nnlsHALSupdt(M,M(:,K),[],100);  
            resultsree(algoix,i) = resultsree(algoix,i) + norm(M-M(:,K)*V,'fro')/nM*100/iter;  
        end
    end
    fprintf('%1.0f...',i);
    if mod(i,10) == 0, fprintf('\n'); end
end

% Display performances
scriptfig; 
% Running time
disp('***********************************');
disp('Total running time: ');  
fprintf('SPA  = %2.2f, \nSNPA = %2.2f, \nXRAY = %2.2f. \n', timings(1), timings(2), timings(3))
disp('***********************************');
% Robustness
disp('Robustness: ');  
for i = 1 : nalgo
    res = results(i,:) >= 1-1e-16; 
    if res == 1
        rob(i) = delta(end);
    else
        [a,b] = min(res);  
        if b == 1, rob(i) = 0; else rob(i) = delta(b-1); end
    end
end
fprintf('SPA  = %1.1d, \nSNPA = %1.1d, \nXRAY = %1.1d. \n', rob(1), rob(2), rob(3)); 
disp('***********************************'); 