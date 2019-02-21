% Testing Lasso HALS for sparse NMF
% Question of reviewer 2 about Figure 2

% Let us try to check whether there is indeed a problem. Let's try
% - giving the right dictionary as input (but wrong coefficients) for sanity check
% - varying the lambda in the paper's setting

% We study the most favorable setting where samples have a dimension d=125 much
% larger than the sparsity level k=3; However the rank is still r=4.

% In addition to studying error on dictionary, contrary to the paper, we
% also investigate the support recovery (i.e. zeros in B).

clc
clear variables
%close all

% ensure current directory is in Matlab PATH

% THE FOLLOWING CODES ARE CURTOSY OF NICOLAS GILLIS:
% - nnlsHALSupdt, nnlsHALSupdtv2
% - HALSacc

% Partial Reproducibility 
%rng(0)

d=10;

%% Parameters
dim = [d 200];
r = 4;
k = 3;% choose 2 or 3 
Nm = 10;

% Storing results
L=8; % Number of lambda tested
stock = zeros(Nm,L);
stockQ = zeros(Nm,L);

% MP feasability
mu = zeros(Nm,1);

%% Data generation


for n=1:Nm

% implementation detail
dim(2) = 200;
    
% random dictionary A 
A = rand(dim(1),r); 
A = col_norm(A); % for dictionary methods
% Random coefficients
B = rand(r,dim(2));
% Project B by block on hyperplanes
for l=0:r-k-1 % zero delta_x position loop
    for p=1:r % initial x zero position loop
        for i=1:(dim(2)/r) % y loop
            B(mod(p+l-1,r)+1,i+dim(2)*(p-1)/r)=0;
        end
    end
end
    
% adding zeros on other segments for uniqueness if k=2
if k==2
   addB = rand(r,50);
   addB(2,:) = 0;
   addB(4,:) = 0;
   B = [B,addB];
   addB2 = rand(r,50);
   addB2(1,:) = 0;
   addB2(3,:) = 0;
   B = [B,addB2];
   
   dim(2) = 300;
end

M = A*B;


%% Running the decompositions

% Initialization
    % Random
%A0 = rand(dim(1),r); A0 = col_norm(abs(A0));
B0 = rand(r,dim(2));

% Oracle for sanity check
A0 = A;
%B0 = B;

% Scaling the intialization
% Scaling, p. 72 of the thesis of N.Gillis
term1 = M*B0'; term2 = B0*B0'; 
scaling = sum(sum(term1.*A0))/sum(sum( term2.*(A0'*A0) )); B0 = B0*scaling; 


%% NMF lasso

% Plain lasso nmf with HALS 
cnt = 1;
perfs  = zeros(r,L); % second arg = number of methods
perfsQ = zeros(L,1);
for lambda = [1 0.5 0.2 1e-1 1e-2 2e-3 1e-3 1e-4]
%for lambda = [ 1e-1 2e-2 1e-2 2e-3 1e-3 5e-4 1e-4]
iter_max_exact = 10000;

% Case 1 : initial A0 is the true A.
% The 'review' version has no acceleration, no normalization and starts by updating B
% (instead of A).
%[Ashals, Bshals, err_shals] = sHALS_review(M, A0, lambda, B0, iter_max_exact); 
%stock_err{n,cnt}=err_shals;
% Case 2 : initial A0 is random (cf ICASSP paper)
%[Ashals, Bshals, err_shals] = sHALSacc(M, A0, lambda, B0, 0.5, 0.01, iter_max_exact);
% One can observe indeed a drop in performances

% First step: fixed (known) dictionary, we only perform sparse coding
UtU = A'*A; UtM = A'*M;
Bshals = HALSupdt_sp_review(B0,lambda,UtU,UtM,iter_max_exact);
Ashals = A0;
err_shals = norm(M-Ashals*Bshals,'fro')^2/norm(M,'fro')^2;
stock_err(n,cnt)=err_shals;
%% Comparing atom estimation performances

% Normalization of the outputs
Ashals = col_norm(Ashals);

% Optimal permutation
perm_shals = munkres(-A'*Ashals); Ashals = Ashals(:,perm_shals);

% Computing error, counting wrong outputs
%sHALS
perfsQ(cnt) = norm(Ashals - A, 'fro')/sqrt(r);

for p=1:r 
    perfs(p,cnt) = (norm(Ashals(:,p) - A(:,p))< 10^(-3)); % Fairer since already very close
    %perfs(p,cnt) = (norm(Ashals(:,p) - A(:,p))< 10^(-4));
end
cnt = cnt+1;
end
stock(n,:) = sum(perfs,1)/r;
stockQ(n,:) = perfsQ';
end


% print performances qualitative (to be averaged over N)
fprintf('\n Atom estimation performance of various lambdas, in %% \n')

display(sum(stock,1)/Nm*100)
rec = sum(stock,1)/Nm*100;
recQ = sum(stockQ,1)/Nm*100;