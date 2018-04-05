%------------------------------D_urban--------------------------------%
% This file tests the Dictionnary ANLS algorithm with hyperspectral data URBAN.
%
%
% List of updates                 -     31/10/2016  -     J. E. Cohen 
%                                       Creation of the file
%                                 -     05/07/2017  -     J. E. Cohen
%                                       Multi-Dictionary version
%                                 -     24/10/2017  -     J. E. Cohen
%                                       Added H2NMF+M2PALS
%-------------------------------------------------------------------------%

%--------------Cleaning the workspace and loading functions---------------%
close     all;
clear     variables;
%-------------------------------------------------------------------------%

% Input your data, second line is the dictionnary
load('Urban.mat')
T         =     A'; %Remark : Do not normalize DATA !!
% Self-dictionnary
D         =     T;
DATA      =     reshape(T,162,307,307);

% Choose dictionary generation method
meth = 'Manual'; % Manual, Unsup, Semisup

% Dictionaries selection

switch(meth)

    case('Manual')
    % Manual
k       = 4;
load('Dm_Urban')
%[ Dm,d,pos ] = D_select( DATA,k );
%[ Dm,d,pos ] = D_select( DATA,k ,indice_mapped);
    
    case('Unsup')

    % Automatic
        % Unsupervised H2NMF
    R_H2NMF = 6;
    d = ones(1,R_H2NMF);
    [idK, ~, ~] = hierclust2nmf(T,R_H2NMF); 
    Dm = cell(1,R_H2NMF);
    corr_table = cell(1,R_H2NMF);
for i=1:R_H2NMF
   indexes = 1:94249;
   index_r = indexes(idK==i);
   Dm{i} = T(:,index_r);
   % Correspondance table
   corr_table{i} = [1:size(index_r,2);index_r]; 
end
            %Printing the segmentation
        label_im = label2rgb(reshape(idK,307,307));
        figure, imagesc(label_im)
        
    case('Semisup')
        % Semi-supervised Segmentation + M2PALS 
        % credits:   Li, J.; Bioucas-Dias, J. M.; Plaza, A.; , "Spectral–Spatial Hyperspectral
        % Image Segmentation Using Subspace Multinomial Logistic Regression and 
        %  Markov Random Fields," Geoscience and Remote Sensing, IEEE Transactions 
        %  on , vol.PP, no.99, pp.1-15, 0  doi: 10.1109/TGRS.2011.2162649
        %
        % To use, download codes at
        % http://www.umbc.edu/rssipl/people/aplaza/ and the graph-cut c++
        % code at http://www.wisdom.weizmann.ac.il/~bagon/matlab.html
        % Graph cut has to be compiled with mex, provided a c++ compiler is
        % configured (GCC for instance). 
      
    % Training data determination
            % hand selection of labeled zones
        k=6;
        load('Dm_Urban3')
        %[ Dm,d,pos ] = D_select( DATA,6 );
            
        
            % Denoising
        [ut,st] = eig(T*T');
        T_pca = ut(:,1:160)'*T;
            % Recovering the pixel indexes
        train = [];
        for i=1:k
            % Conversion from rect position to index
            index1 = (pos(i,1)-1)*307 + pos(i,2);
            index_train = index1;
            for l=1:pos(i,3)
                index_train = [index_train,(index1+(l-1)*307)+(1:pos(i,4))];
            end
            % Associating each pixel in each zone with a label
            train = [train,[index_train;i*ones(1,size(index_train,2))]];
        end
          train_dat = T_pca(:,train(2,:));           
          [v_seg,d_seg]   = eig(train_dat*train_dat'/size(train,2));
          d_seg = diag(d_seg);
          dtrue = d_seg(25); %68
          g = [];
          g_all= [];
            
            % Using segmentation as in exemple file
        for k_iter = 1:k
            train_k = train(2,:) == k_iter;
            train_dat_k = train_dat(:,train_k);     
            n_k   = size(train_k,2);
            [v_seg,d_seg]   = eig(train_dat_k*train_dat_k'/n_k);
            d_seg = diag(d_seg);
            
            sub_sp = d_seg<dtrue;
            sub(k_iter) = sum(sub_sp);
            tau = sub(k_iter);
            
            P = v_seg(:,1:tau)*v_seg(:,1:tau)'*train_dat;
            n_train = size(train,2);
            gall = zeros(1,n_train);
            
            for num_k = 1:n_train
                gall(num_k) = sqrt(P(:,num_k)'*P(:,num_k));
            end
            g = [g; gall];            
            n_all = 307*307;
            ggall=zeros(1,n_all);
            P_all = v_seg(:,1:tau)*v_seg(:,1:tau)'*T_pca;
            
            for iter_all = 1:n_all
                ggall(iter_all) = sqrt(P_all(:,iter_all)'*P_all(:,iter_all));
            end
            g_all = [g_all;ggall];
        end
          
            % regularization parameter
        lambda = eps;
        w = subspace_classifier(-g,train(2,:),lambda);

            % compute the probablity
        p= subspace_mlogistic(w,-g_all);
       [maxp,class] = max(p);
       
            % Segmentation by graph_cut
        Dc = reshape((log(p+eps))',[307 307 k]);
        Sc = ones(k) - eye(k);        
        beta = 2; % increase for regularity (4 ?)
        gch = GraphCut('open', -Dc, beta*Sc);
        [gch seg] = GraphCut('expand',gch);
        gch = GraphCut('close', gch);
        class = double(seg(:));
        
            % Printing the segmentation
        label_im = label2rgb(reshape(class,307,307));
        figure, imagesc(label_im)
        
            % Storing into Dm
        Dm = cell(1,k);
        corr_table = cell(1,k);
        for i=1:k
        indexes = 1:94249;
        index_r = indexes(class==(i-1));
        Dm{i} = T(:,index_r);
        % Correspondance table
        corr_table{i} = [1:size(index_r,2);index_r]; 
        end    
            % 2 components can be taken in each zone
        d = 2*ones(1,k);
        
end
%% -----------------------------Model parameters----------------------------%
% 3-way block dimensions
[m,n]     =     size(T);
% Number of components
R         =     6;
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
keepy = mod(K_omp-1,307)+1;
keepx = floor((K_omp-1)/307)+1;
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
indice_2D = mod(K_omp_true-1,307)+1;
indice_1D = floor((K_omp_true-1)/307)+1;
indice_mapped = [indice_1D;indice_2D];

    case('Semisup')
% With Semi-supervised Automated dictionary
for i=1:R
    K_omp_true(i)  = corr_table{K_omp(1,i)}(2,K_omp(2,i));
end
indice_2D = mod(K_omp_true-1,307)+1;
indice_1D = floor((K_omp_true-1)/307)+1;
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
    imagesc(reshape(B_omp_old(:,i),307,307))%307 x n
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
    subplot(2,R,i)
    imagesc(reshape(B_omp(:,i),307,307))%307 x n
    xlabel(norms(i))
    subplot(2,R,R+i)
    %plot(A_omp(:,i))
    %hold on
    plot(A_omp(:,i),'r')
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
h=reshape(sum(T,1),307,307);
imagesc(h)
hold on
for i=1:k
    rectangle('Position',pos(i,:),'EdgeColor','k','LineWidth',1.5)
end
elseif strcmp(meth,'Semisup')
% for semi-supervised segmented   
imagesc(label_im)
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
hold on
plot(keepx,keepy,'xblack','MarkerSize',20,'LineWidth',1.5)
title('Residuals and selected pixels')

%----------- Spectral Variation plot ----------
figure
h2=reshape(abs(sum(B_omp,2)-1),307,307);
imagesc(h2)
title('Variation around sum to 1')
