%------------------------------D_SanDiego------------------------------%
% This file tests the Dictionnary ANLS algorithm with hyperspectral data SanDiego.
%
% External functions:
%
% - kr.m              : Khatri-Rao product;
% - D_anls2.m         : OMP-ANLS algorithm;
% - D_anls_flex.m     : Flex-OMP-ANLS algorithm;
% - D_admm.m          : ADMM algorithm;
% - amb_correct.m     : Matching Permutation for two CP models.
%
% List of updates                 -     27/01/2017  -     J. E. Cohen 
%                                       Creation of the file
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
clc;
clear     all;
close     all;
%-------------------------------------------------------------------------%

% Input your data, second line is the dictionnary
load('SanDiego.mat')

coln      =     400;
T         =     max(A((1+400*0):(400*coln),:)',0);

% Self-dictionnary
D         =     T;
%-----------------------------Model parameters----------------------------%
% 3-way block dimensions
[m,n]     =     size(T);
[md,nd]   =     size(D);
% Number of components
R         =     10;
% Number of iterations
iter      =     50;


%-------------------------------Initialization----------------------------------% 
r = 2;
tic

if r == 1
        % Initialize with SPA
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
        % FGNSR-100 
        clus_num = 100; toc;
        [IDX,~,K0] = hierclust2nmf(T,clus_num); 
        DD       =  T(:,K0);
        for ki = 1 : max(IDX)
            w(ki) = sqrt(sum(IDX==ki));
        end
        tic
        [~,K] = fgnsr(DD*diag(w),R,ones(size(DD,2),1),'maxiter',1000,'verbose', false, 'proj', 'parproj', 'col_weights', ones(size(DD,2),1)); 
        disp(' FGNSR-100 Error ')
        K     = K0(K);  
end

A_init    =     D(:,K);
B_init    =     nnlsHALSupdt(T,A_init,[],500);   
toc
err_init  =     sqrt(sqFnorm(T-A_init*B_init))/sqrt(sqFnorm(T))
B_init    =     nnlsHALSupdt(T,A_init,[],1);   % NNLS: M approx M(:,K) * Vspa
B_init    =     B_init';


%--------------------------Dictionnary decomposition----------------------%
opts{2} = 'true';opts{3}={'NMFinit'};

% MP Glouton 
tic
[A_omp,B_omp,K_omp,err_omp]=D_nmf_OMP(T,iter,B_init,K,D,opts);
toc
% MP Glouton Nicolas
%[A_omp,B_omp,K_omp,err_omp]=D_nmf_SMP(T,iter,B_init,K,D,opts);

% MP Flex Glouton
tic
[A_omp,B_omp,K_omp,err_omp]=D_nmf_FlexMP(T,250,B_init,K,D,100,0.1,opts);

% Fast Gradient 
%[B_omp,A_omp,K_omp,err_grad] = D_nmf_fastgrad(T',iter,B_init,A_init,D,opts);

toc

Mres = abs(T-A_omp*B_omp');

indice_2D = mod(K_omp,400);
indice_1D = floor(K_omp/400);
indice_mapped = [indice_1D;indice_2D];

norms   =   sum(B_omp.^2);
norms   =   floor(1000*norms/sum(norms))/10;
B_omp   =   col_norm(B_omp);
figure
for i=1:R
    subplot(2,R,i)
    imagesc(reshape(B_omp(:,i),400,coln))%400 x n
    xlabel(norms(i))
    subplot(2,R,R+i)
    plot(A_omp(:,i))
    hold on
    plot(D(:,K_omp(i)),'r')
end

% normsi  =   sum(B_init.^2);
% normsi  =   floor(1000*normsi/sum(normsi))/10;
% %B_init  =   col_norm(B_init);
% figure
% for i=1:R
%     subplot(2,R,i)
%     imagesc(reshape(B_init(:,i),400,coln))%400 x n
%     xlabel(normsi(i))
%     subplot(2,R,R+i)
%     plot(A_init(:,i))
% end
% 
figure
imagesc(reshape(Mres(50,:),400,coln))
hold on
plot(indice_1D,indice_2D,'+red')