% Comparison test file for K-sparse MF (sparsity on B)
% Final ICASSP version

clc
clear variables
%close all
% ensure current directory is in Matlab PATH

% WARNING: THIS CODE USES EXTERNAL CODE:
% - ksubspace : https://fr.mathworks.com/matlabcentral/fileexchange/37353-k-subspaces?focused=5236203&tab=function
% - sNNOMP: ask the authors of 
% `` An optimized version of non-negative OMP'', GRETSI 2017

% ALSO, THE FOLLOWING CODES ARE CURTOSY OF NICOLAS GILLIS:
% - nnlsHALSupdt, nnlsHALSupdtv2
% - HALSacc

% The following code is provided without explicit author consent
% - nnlsm_activeset + solveNormalEqComb


% Partial Reproducibility 
rng(0)

% Main loop
dlist = [4,5,10,20,25,50,125];
counter = 0;
%lambda = [0.1 0.1 0.1 0.01 0.01 0.01 0.01]; % for k=3, has been chosen so that LASSO HALS works best
lambda = [0.1 0.2 0.1 0.1 0.1 0.1 0.1]; % for k=2
for d=dlist
counter = counter+1;

%% Parameters
%constraints = 'NN'; % nothing, or 'NN'
noise = 0;
dim = [d 200];
r = 4;
k = 2;% choose 2 or 3 
Nm = 2;
iter_max_exact = 1000;

% Storing results
L = 7;
stock{counter} = zeros(Nm,L);
stockQ{counter} = zeros(Nm,L);

% MP feasability
mu = zeros(Nm,1);

%% Data generation


for n=1:Nm

% implementation detail
dim(2) = 200;
    
% random A 
A = rand(dim(1),r); 
A = col_norm(A); % for dictionary methods

% Project B by block on hyperplanes
B = rand(r,dim(2));
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
%Mtrue = A*B;
%M = A*B+noise*randn(dim);
 %SNR = 10*log10(norm(Mtrue,'fro')^2/norm(M-Mtrue,'fro')^2);

%% Theoretical guarenties at optimum
% 
% % Condition 1: Coherence
% Gram = A'*A;
% Gram(eye(r)==1) = 0;
% mu(n) = max(max(abs(Gram)));
% kmin = (1+1/mu(n))/2;
% if k<kmin
%     fprintf('Coherence condition satisfied at optimum, MP methods may succeed')
% else
%     fprintf('Coherence condition not satisfied, MP methods will fail') 
% end
% 
% % Condition 2: ERC (CNS)
% % TODO, combinatorial

%% Running the decompositions



% Initialization
    % Random
A0 = rand(dim(1),r); A0 = col_norm(abs(A0));
B0 = rand(r,dim(2));
%A0 = [eye(dim(1)),rand(dim(1),r-dim(1))] + 0.1*rand(dim(1),r);
%A0 = [eye(dim(1)),rand(dim(1),r-dim(1))];

% Oracle for sanity check
%A0 = A;
%B0 = B;

    % Pure pixel selection
% SPA init + nnls sur B
%K0 = FastSepNMF(M, r);
%A0 = M(:, K0); A0 = col_norm(A0);
%B0 = nnlsHALSupdt(M, A0, [], 500);

% Scaling the intialization
% Scaling, p. 72 of the thesis of N.Gillis
term1 = M*B0'; term2 = B0*B0'; 
scaling = sum(sum(term1.*A0))/sum(sum( term2.*(A0'*A0) )); B0 = B0*scaling; 

%% ESNA

% Combinatorial K-NMF
    % Decomposition
 %   tic
[Aknmf,Bknmf,errn] = k_s_nmf(M,k,r,iter_max_exact,A0,B0);
 %   toto=toc;
%% NMF lasso

% NMF only
[Anmf, Bnmf, err_nmf] = HALSacc(M, A0, B0, 0.5, 0.01, iter_max_exact);

% Lasso nmf with HALS 
%lambda = 1e-3; % old value
[Ashals, Bshals, err_shals] = sHALS_review(M, A0, lambda(counter), B0, 10*iter_max_exact);

%% NNOMP

% Truncated active set 
[Atas, Btas, err_tas] = as_k_s_nmf(M,k,iter_max_exact,A0,B0);

% ---------------- COMMENTED ----------
% % Alternated SNNOMP (credits to Thi Than Nguyen et. al.)
% [Asnn, Bsnn, err_snn] = snnomp_k_s_nmf(M,k,500,A0,B0);
% -------------------------------------

%% Subspace clustering

% K-summits (aka NOLRAK) (proposed)
tic
[Aksum, Bksum, affect_ksum,err_ksum] = k_summits(M,k,iter_max_exact,A0,B0);
titi=toc;

% ----------------- COMMENTED --------------------
% K-planes (credits to Masayuki Tanaka) 

% [Um,Sm,Vm] = svd(M); Um = Um(:,1:r);
% Mproj = Um'*M;
% A0proj = Um'*A0;
%     %Initial SS
%     SS0 = zeros(r-1,r,r);
%     SS00 = cell(r,1);
%     for i=1:r
%         SS0(:,:,i) = A0proj(:,setdiff(1:r,i))';
%         [U,S,V] = svd(A0proj(:,setdiff(1:r,i)),'econ');
%         SS00{i} = U';
%     end
%     % Clustering
% [IDX,SS] = ksubspaces(Mproj',r,r-1,SS0,500);
% IDX_KH = IDX';
% 
% % Computing the estimated atoms by subspace intersection
%     % 1. normals to the subspaces
%     N = zeros(r,r);
%     for p=1:r
%         np = null(SS(:,:,p));
%         N(:,p) = np;
%     end
%     
%     % 2. intersecting the hyperplanes
%     N = N';
%     Akp = zeros(r,r);
%     for p=1:r
%        ap = null(N(setdiff(1:r,p),:));
%        Akp(:,p) = ap;
%     end
%     
% Akp = Um*Akp; 
%------------------------------------------------------------
    
%% Comparing atom estimation performances

% Normalization of the outputs
Aknmf = col_norm(Aknmf);
Anmf = col_norm(Anmf);
%Akp = col_norm(abs(Akp));
Ashals = col_norm(Ashals);
Atas = col_norm(Atas);
%Asnn = col_norm(Asnn);
Aksum = col_norm(Aksum);


% Optimal permutation
perm_nmf = munkres(-A'*Anmf); Anmf = Anmf(:, perm_nmf);
perm_en = munkres(-A'*Aknmf); Aknmf = Aknmf(:,perm_en);
%perm_ekp = munkres(-A'*Akp); Akp = Akp(:,perm_ekp);
perm_shals = munkres(-A'*Ashals); Ashals = Ashals(:,perm_shals);
perm_tas = munkres(-A'*Atas); Atas = Atas(:,perm_tas);
%perm_snn = munkres(-A'*Asnn); Asnn = Asnn(:,perm_snn);
perm_ksum = munkres(-A'*Aksum); Aksum = Aksum(:,perm_ksum);

% Computing error, counting wrong outputs
% 1- KNMF [proposed]
% 2- Kplanes (commented out)
% 3- NMF
% 4- sHALS
% 5- truncated activetset NNOMP 
% 6- support NNOMP (commented out)
% 7- Ksum [proposed]

perfs  = zeros(r,L); % second arg = number of methods
perfsQ(1) = norm(Aknmf - A, 'fro')/sqrt(r);   
%perfsQ(2) = norm(Akp - A, 'fro')/sqrt(r);  
perfsQ(3) = norm(Anmf - A, 'fro')/sqrt(r);
perfsQ(4) = norm(Ashals - A, 'fro')/sqrt(r);
perfsQ(5) = norm(Atas - A, 'fro')/sqrt(r);
%perfsQ(6) = norm(Asnn - A, 'fro')/sqrt(r);
perfsQ(7) = norm(Aksum - A, 'fro')/sqrt(r);

prec = 10^(-4);
for p=1:r 

    perfs(p,1) = (norm(Aknmf(:,p) - A(:,p))< prec);   
%    perfs(p,2) = (norm(Akp(:,p) - A(:,p))< prec);
    perfs(p,3) = (norm(Anmf(:,p) - A(:,p))< prec);
    perfs(p,4) = (norm(Ashals(:,p) - A(:,p))< prec);
    perfs(p,5) = (norm(Atas(:,p) - A(:,p))< prec);
%    perfs(p,6) = (norm(Asnn(:,p) - A(:,p))< prec);
    perfs(p,7) = (norm(Aksum(:,p) - A(:,p))< prec);
end

stock{counter}(n,:) = sum(perfs,1)/r;
stockQ{counter}(n,:) = perfsQ;
end

end

% printing out results

%fprintf('ENSA, KP, NMF, sHALS, truncated active-set NNOMP, SuNNOMP, NOLRAK')
fprintf('ENSA, NMF, sHALS, truncated active-set NNOMP, NOLRAK')

% print performances qualitative (to be averaged over N)
fprintf('\n Atom estimation performance of various methods, in %% \n')
for counter=1:size(dlist,2)
display(sum(stock{counter},1)/Nm*100)
plotcurves(counter,:) = sum(stock{counter},1)/Nm*100;
plotcurvesQ(counter,:) = sum(stockQ{counter},1)/Nm*100;
end

% print performance quantitative WHEN FAILURE ONLY (when no mistake, 0%)
%fprintf('\n Relative estimation performance on A of various methods, failure cases only, in %% \n')
%for counter=1:6
%display(sum(stockQ{counter},1)./(Nm*ones(1,L)-sum(stock{counter},1)))
    %%disp(std(stockQ{counter},1))
%end

% Plots
namearray={'LineStyle', 'Marker', 'Color'};
valuearray={'-','d','k';...  
    '--','o','k';...
    '--','x','k';...
    '--','d','r';...
    '-','o','r';...
    '--','x','r';...
    '-','d','b'};
%figure

if k==3
    subplot(221)
else
    subplot(223)
end
fig=semilogx(dlist,plotcurves,'LineWidth',2,'MarkerSize',10)
set(fig, namearray, valuearray);
if k==3
legend('ESNA(proposed)', 'KSub','NMF-HALS', 'Lasso-HALS','a.set NNOMP', 'SuNNOMP', 'NOLRAK(proposed)')
end
xlabel('d')
ylabel('Quasi perfect reconstruction of D (%)')
xticks(dlist)


if k==3
    subplot(222)
else
    subplot(224)
end
fig=semilogx(dlist,plotcurvesQ,'LineWidth',2,'MarkerSize',10)
set(fig, namearray, valuearray);
%legend('ESNA(proposed)', 'KSub','NMF-HALS', 'Lasso-HALS','a.set NNOMP', 'SuNNOMP', 'NOLRAK(proposed)')
xlabel('d')
ylabel('relative MSE on D')
xticks(dlist)