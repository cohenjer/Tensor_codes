%------------------------------D_cuprite--------------------------------%
% This file tests the Dictionnary ANLS algorithm with hyperspectral data Cuprite.
%
%
% List of updates                 -     07/12/2017  -     J. E. Cohen 
%                                       Creation of the file
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
close     all;
clear     variables;
%-------------------------------------------------------------------------%

% Input your data, second line is the dictionnary
%load('cuprite_ref.mat')
T         =     max(x,0); %Remark : Do not normalize DATA !!
clear x;
% Self-dictionnary
D         =     T;
DATA      =     reshape(T,188,250,191);

% Choose dictionary generation method
meth = 'Unsup'; % Manual, Unsup

% Dictionaries selection

switch(meth)

    case('Manual')
    % Manual
k       = 8;
%load('Dm_Urban')
[ Dm,d,pos ] = D_select( DATA,k );
%[ Dm,d,pos ] = D_select( DATA,k ,indice_mapped);
    
    case('Unsup')

    % Automatic
        % Unsupervised H2NMF
    R_H2NMF = 12;
    d = ones(1,R_H2NMF);
    [idK, ~, ~] = hierclust2nmf(T,R_H2NMF); 
    Dm = cell(1,R_H2NMF);
    corr_table = cell(1,R_H2NMF);
for i=1:R_H2NMF
   indexes = 1:47750;
   index_r = indexes(idK==i);
   Dm{i} = T(:,index_r);
   % Correspondance table
   corr_table{i} = [1:size(index_r,2);index_r]; 
end
            %Printing the segmentation
        label_im = label2rgb(reshape(idK,250,191));
        figure, imagesc(label_im)

        
end
%% -----------------------------Model parameters----------------------------%
% 3-way block dimensions
[m,n]     =     size(T);
% Number of components
R         =     12;
%R         =     sum(d);
% Number of atoms
K         =     size(n,2);
% Number of iterations
iter      =     50;
            

%-------------------------------Initialization----------------------------------% 
% To be modified
r = 4; %4 for Unsup, Semisup
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
        w        =  zeros(max(IDX),1);
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
err_init  =     sqrt(sqFnorm(T-A_init*B_init))/sqrt(sqFnorm(T));
disp(err_init)
B_init    =     nnlsHALSupdt(T,A_init,[],1);   % NNLS: M approx M(:,K) * Vspa
B_init    =     B_init';

%--------------------------Dictionnary decomposition----------------------%

% OPTIONS
opts{2} = 'true';opts{3}='NMFinit';opts{4}='eucl';

% MP ALS 
tic
[A_omp,B_omp,K_omp,err_omp]=D_nmf_OMP(T,iter,B_init,K,D,opts);
toc
% Getting indexes with normal dictionary
keepy = mod(K_omp-1,250)+1;
keepx = floor((K_omp-1)/250)+1;
K_omp_old = K_omp;
keep  = [keepx;keepy]; % Indexes of atoms selected by MPALS
B_omp_old = B_omp;
A_omp_old = A_omp;

% MP ALS Multi Dico
% % --- test zone ----
% for i=1:R
%     %Dm{i} = A_omp_old(:,R-i+1);
%     Dm{i} = [T(:,randperm(94249,45)),A_omp_old(:,R-i+1),T(:,randperm(94249,50))];
% end
% d = [1 1 1 1 1 1];
%opts{3} = 'NMFinit';

opts{1} = 'NN';opts{2}='NMFinit';opts{3}='eucl';

tic
[A_omp,B_omp,K_omp,err_mult_omp]=Dmult_nmf_OMP(T,100,B_init,A_init,Dm,d,opts);
toc
Mres = sum(abs(T-A_omp*B_omp'),1);

%------------- Getting the position of extracted pure spectra -------

% With multiple dictionaries
switch(meth)
    
    case('Manual')
indice_1D = zeros(1,R);
indice_2D = zeros(1,R);
for i=1:R
   ymax = pos(K_omp(1,i),4);
   indice_yi = mod(K_omp(2,i)-1,ymax+1); indice_xi = floor((K_omp(2,i)-1)/(ymax+1));
   indice_1D(i) = pos(K_omp(1,i),1)+indice_xi;
   indice_2D(i) = pos(K_omp(1,i),2)+indice_yi;
end
indice_mapped = [indice_1D;indice_2D];


    case('Unsup')
% With H2NMF Automated dictionary
for i=1:R
    K_omp_true(i)  = corr_table{K_omp(1,i)}(2,K_omp(2,i));
end
indice_2D = mod(K_omp_true-1,250)+1;
indice_1D = floor((K_omp_true-1)/250)+1;
indice_mapped = [indice_1D;indice_2D];

end

% ------------ Plotting Results and Processing ------------

%--------- OMPALS plots ----------- 
norms   =   sum(B_omp_old.^2); %% CORR : extraire la norme de A
norms   =   floor(1000*norms/sum(norms))/10;
%B_omp_old   =   col_norm(B_omp_old);
figure 
for i=1:R
    subplot(2,R,i)
    imagesc(reshape(B_omp_old(:,i),250,191))%307 x n
    xlabel(norms(i))
    subplot(2,R,R+i)
    plot(A_omp_old(:,i),'r')
end

%--------- OMP ALS Mult plots ----------- 
norms   =   sum(B_omp.^2);
norms   =   floor(1000*norms/sum(norms))/10;
%B_omp   =   col_norm(B_omp);
figure 
for i=1:R
    subplot(2,R/2,i)
    imagesc(reshape(B_omp(:,i),250,191))%307 x n
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    %xlabel(norms(i))
    %subplot(2,R,R+i)
    %plot(A_omp(:,i))
    %hold on
    %plot(A_omp(:,i),'r')
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
end

%--------- Init plots ----------- 
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
%

%------------ Residual Plots ------------
figure
if strcmp(meth,'Manual')
% h=reshape(Mres,307,307); residuals
h=reshape(sum(T,1),250,191);
imagesc(h)
hold on
for i=1:k
    rectangle('Position',pos(i,:),'EdgeColor','k','LineWidth',1.5)
end
else
% for segmented images
imagesc(label_im) % obtained at the beginning of the test
end

hold on
plot(indice_1D,indice_2D,'+black','MarkerSize',20,'LineWidth',1.5)
%hold on
%plot(keepx,keepy,'xblack','MarkerSize',20,'LineWidth',1.5)
%title('Residuals and selected pixels')

%----------- Spectral Variation plot ----------
figure
h2=reshape(abs(sum(B_omp,2)-1),250,191);
imagesc(h2)
title('Variation around sum to 1')
