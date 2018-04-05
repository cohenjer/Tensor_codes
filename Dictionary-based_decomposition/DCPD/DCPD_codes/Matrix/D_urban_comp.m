%------------------------------D_urban--------------------------------%
% This file tests the Dictionnary ANLS algorithm with hyperspectral data URBAN.
%
% External functions:
%
% - kr.m              : Khatri-Rao product;
% - D_anls2.m         : OMP-ANLS algorithm;
% - D_anls_flex.m     : Flex-OMP-ANLS algorithm;
% - D_admm.m          : ADMM algorithm;
% - amb_correct.m     : Matching Permutation for two CP models.
%
% List of updates                 -     31/10/2016  -     J. E. Cohen 
%                                       Creation of the file
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
clear     variables;
close     all;
%-------------------------------------------------------------------------%

% Input your data, second line is the dictionnary
load('Urban.mat')

% Normalize image
A         =     A/max(A(:));

% Extracting the top left part of the image
siz       =     150;
Im        =     reshape(A,307,307,162);
Im        =     Im(1:siz,1:siz,:);
T         =     reshape(Im,siz^2,162)';
%T         =     A((1+307*0):(307*100),:)'; %size has been reduced for computing

% Self-dictionnary
%D         =     T(:,randperm(siz*siz,200));
% Number of atoms
%K         =     size(D,2);

% Presentation plots
figure
h=imagesc(reshape(T(10,:),siz,siz));

%-----------------------------Model parameters----------------------------%
% 3-way block dimensions
[m,n]     =     size(T);
% Number of components
R         =     6;
% Number of iterations
iter      =     20;
            

%-------------------------------Initialization----------------------------------% 
r = 6;
tic

if r == 1 % For now, bugged for r\=5
        K = FastSepNMF(T,R); 
        disp(' SPA Error ')
    elseif r == 2
        % Initialize with VCA
        K = VCA(T,'Endmembers',R); 
        disp(' VCA Error ')
    elseif r == 3
        % Initialize with SNPA
        K = SNPA(T,R);  
        disp(' SNPA Error ')
    elseif r == 4    
        % Initialize with H2NMF
        [~, ~, K] = hierclust2nmf(T,R); 
        disp(' H2NMF Error ')
    elseif r == 5 
        % Random init 
        n = size(D,2); 
        K = randperm(n,R); 
        disp(' Random init Error ')
    elseif r == 6 
        % FGNSR-200 
        clus_num = 200; toc;
        [IDX,~,K0] = hierclust2nmf(T,clus_num); 
        %DD       =  T(:,K0);
        D        =  T(:,K0);
        w        =  zeros(max(IDX),1);
        for ki = 1 : max(IDX)
            w(ki) = sqrt(sum(IDX==ki));
        end
        tic
        [~,K] = fgnsr(D*diag(w),R,ones(size(D,2),1),'maxiter',1000,'verbose', false, 'proj', 'parproj', 'col_weights', ones(size(D,2),1)); 
        disp(' FGNSR-200 Error ')
        %K     = K0(K);  
end

A_init    =     D(:,K);
B_init    =     nnlsHALSupdt(T,A_init,[],500);   
toc
res_init  =     norm(T-A_init*B_init,'fro')/norm(T,'fro');
disp(res_init)
B_init    =     nnlsHALSupdt(T,A_init,[],1);   % NNLS: M approx M(:,K) * Vspa
B_init    =     B_init';


%--------------------------Dictionnary decomposition----------------------%
opts{2} = 'true';opts{3}={'NMFinit'};%opts{3}={''};

% MP Glouton 
[A_omp,B_omp,K_omp,err_omp]=D_nmf_OMP(T,iter,B_init,K,D,opts);
res_omp = norm(T-A_omp*B_omp','fro')/norm(T,'fro');

% GLUP
rho = 10; mu = 1; epsabs = 1E-2; epsrel = 1E-2; % to be optimized
[X_glup,iter_glup]=GLUP(D,D,rho,mu,epsabs,epsrel);
meanX = mean(X_glup,2);
[NX_sort, idx] = sort(meanX,'descend');
%tol = 0.03;
%R_tol = sum(meanX>tol);
R_tol  = R;
idx = idx(1:R_tol);
A_glup = T(:,idx);
B_glup = nnlsHALSupdt(T,A_glup,[],500);
res_glup = norm(T-A_glup*B_glup,'fro')/norm(T,'fro');

% SDSOMP
%[A_sdsomp,K_sdsomp] = SDSOMP(T,R);
%B_sdsomp = nnlsHALSupdt(T,A_sdsomp,[],500); 
[A_sdsomp,B_sdsomp] = SOMP(T,D,R);
res_sdsomp = norm(T-A_sdsomp*B_sdsomp,'fro')/norm(T,'fro');

Mres = sum(abs(T-A_omp*B_omp'),1);
%Mres = abs(T-A_omp*B_omp');

indice_2D = mod(K_omp,siz);
indice_1D = floor(K_omp/siz);
indice_mapped = [indice_1D;indice_2D];

norms   =   sum(B_omp.^2);
norms   =   floor(1000*norms/sum(norms))/10;
B_omp   =   col_norm(B_omp);
figure 
for i=1:R
    subplot(2,R,i)
    imagesc(reshape(B_omp(:,i),siz,siz))%307 x n
    xlabel(norms(i))
    subplot(2,R,R+i)
    %plot(A_omp(:,i))
    %hold on
    plot(D(:,K_omp(i)),'r')
end

% 
% normsi  =   sum(B_init.^2);
% normsi  =   floor(1000*normsi/sum(normsi))/10;
% %B_init  =   col_norm(B_init);
% figure
% for i=1:R
%     subplot(2,R,i)
%     imagesc(reshape(B_init(:,i),307,307))%307 x n
%     xlabel(normsi(i))
%     subplot(2,R,R+i)
%     plot(A_init(:,i))
% end

figure
imagesc(reshape(Mres,siz,siz))
figure(1) % for presentation
hold on
plot(indice_1D,indice_2D,'+red')

