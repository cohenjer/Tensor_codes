%------------------- Synthetic_Compar_Urban.m -------------%
% This script is written to compare the performances of M2PALS with other
% sparse-coding algorithms (SPA, SNPA, SDOMP, GLUP, XRAY, FGNSR). These methods
% are contained in the folder ``Tools''. 
% 
% Because the problem at hand is different for M2PALS and requires a
% hand-selection of regions of interest, this file build a synthetic HSI using
% spectra extracted from the URBAN HSI (get link below). A *** distribution is
% used to generate the abundances, and we pick the dictionaries for M2PALS so
% that they indeed contain d_i pure pixels.
%
% The grid parameters are the noise level applied on the synthesized HSI, the
% size of the selected zones and the number of pure pixels per zones. Ideally,
% the results should not deteriorate too much when increasing window size or
% pure pixel count, and M2PALS is hoped to be robust to noise for small windows
% around a pure pixel. 
%
% References:
% Urban data set:  
% SPA:
% SNPA:
% SDOMP:
% GLUP:
% XRAY:
% FGNSR:
%
% Credits:                                              03/08/2017 J.E.Cohen
%                                                       Creation of the file

%---------------------------------
clc
close all
clear variables
%---------------------------------
% Loading the real HSI URBAN
load('Urban.mat')
T = A'; % T is a 162 x 307^2 matrix. 
R = 6; % We only need R spectra.

% Using SPA to extract R good spectra from the data
K = FastSepNMF(T,R); 
A_synth = T(:,K); % A_synth is the synthetic factor A obtained by extracting
% spectra from T, of size 162 x R


%----------------------------------------------
% Main loop, repeating the experiment with different pure pixels positions and
% abundances.

% Number of repetitions of the experiment
N = 100;
% Initialization of the result stocking table
store = zeros(12,5,N); % 12 is the number of methods, including 4 window sizes
% for M2PALS and H2NMF+M2PALS
store_t = zeros(12,5,N); % Same as store but for computation time.

for n=1:N
    clc
    fprintf('ITERATION NUMBER %d\n',n)
    % Generating the abundances B using a sparsity-inducing distribution (Poisson
    % or Dirichlet)
    size_synth = 200; %1000
    % POISSON
    % B_synth = poissrnd(20,size_synth,R); % Poisson distribution with parameter lambda=1
    % induce sparsity.
    % % sumB = sum(B_synth,2); % Sum of columns of B, for adding at least one non-zero value per column of B
    % % for i=1:size_synth 
    % %     if sumB(i) == 0 
    % %         temp_index = randperm(R,2);
    % %         B_synth(i,temp_index) = 0.5; % Adding a non-zero value, for non-empty & non-pure pixels
    % %     end
    % % end
    % % B_synth = B_synth./repmat(sum(B_synth,2),1,R); % normalizing row-wise
    % DIRICHLET
    B_synth = sample_dirichlet(ones(R,1),size_synth); % Dirichlet distribution with
    % parameters alpha_i = alpha. Abundances sum to one and may be sparse.
    % Note that in this data, there is no notion of continuity of
    % similar pixels in the image. This is not a problem for usual sparse coding
    % methods which do not make use of spatial information, but this may bias
    % results for M2PALS since a dictionary for a real HSI will have more
    % correlation than in this synthetic data set.
    % Adding pure pixels
    pp_index = randperm(size_synth,R); % Choosing the pure pixels location
    % randomly. May be modified to place them on a grid, for easier computation and
    % visualization later on.
    B_synth(pp_index,:) = eye(R);

    % Generating the simulated HSI
    M_synth = A_synth*B_synth';

    % ---------------------------
    % Second loop, over noise level

    % Definition of the various values of SNR
    SNR = [0,10,20,30,50];
    count = 0;

    % looping over SNR
    for i=SNR
        count = count+1;
        fprintf('--------- SNR %d -------\n',i)
        var = 10^-(i/10)*norm(M_synth,'fro')^2/R/size_synth; % Computing the
        % noise variance knowing the SNR.
        % Data = M_synth + sqrt(var)*randn(size(M_synth));
        Data = max(M_synth + sqrt(var)*randn(size(M_synth)),0); % without negative values.
        Data = Data/max(max(Data)); % Normalisation of the data, useful of GLUP
        % The SNR is wrong however in this second line

        % Unmixing using sparse coding methods
            % Finding 200 good pixels
    %         clus_num = 200; 
    %         [IDX,~,K0] = hierclust2nmf(Data,clus_num); 
    %         DD       =  Data(:,K0);
            DD = Data;

            % SPA
            tic
            K = FastSepNMF(Data,R); 
            store(1,count,n) = R - length(intersect(K,pp_index));
            %fprintf('SPA %d \n',R - length(intersect(K,pp_index)))
            store_t(1,count,n) = toc;

            % SNPA
            tic
            K = SNPA(Data,R);
            store(3,count,n) = R - length(intersect(K,pp_index));
            %fprintf('SNPA %d \n',R - length(intersect(K,pp_index)))
            store_t(3,count,n) = toc;

            % MPALS (SNPA init) % 
            tic
            B_init = nnlsHALSupdt(Data,Data(:,K),[],500)';
            K_init = K;
            opts{2} = 'true'; opts{3} = 'NMFinit';
            opts{4} = 'eucl';
            [~,~,K] = ...
            D_nmf_OMP(Data,1000,B_init,K_init,Data,opts); 
            store(2,count,n) = R - length(intersect(K,pp_index));
            %fprintf('MPALS %d \n', R - length(intersect(K,pp_index)))
            store_t(2,count,n) = toc;


            % FGNSR-100
            tic
            [~,K] = fgnsr(DD,R,ones(size(DD,2),1),'maxiter',1000,...
            'verbose', false, 'proj', 'parproj', 'col_weights', ones(size(DD,2),1)); 
            store(4,count,n) = R - length(intersect(K,pp_index));
            %fprintf('FGNSR %d \n', R - length(intersect(K,pp_index)))
            store_t(4,count,n) = toc;

            % GLUP
            tic
            rho = 10; mu = 1; epsabs = 1E-2; epsrel = 1E-2; 
            [X_glup,iter_glup]=GLUP(Data,DD,rho,mu,epsabs,epsrel);
            meanX = mean(X_glup,2);
            [NX_sort, idx] = sort(meanX,'descend');
            K = idx(1:R);
            store(5,count,n) = R - length(intersect(K,pp_index));
            %fprintf('Glup %d \n', R - length(intersect(K,pp_index)))
            store_t(5,count,n) = toc;
            
            % N-FINDR
            % PCA will be performed automatically by N-FINDR
            tic
            [~,K] = Nfindr(Data,R);
            store(6,count,n) = R - length(intersect(K,pp_index));    
            %fprintf('NFINDR %d \n', R - length(intersect(K,pp_index)))
            store_t(6,count,n) = toc;

            % SDOMP
            tic
            [~,~,K] = SOMP(Data,Data,R);
            store(7,count,n) = R - length(intersect(K,pp_index));
            %fprintf('SNPA %d \n', R - length(intersect(K,pp_index)))
            store_t(7,count,n) = toc;


            % M2PALS
            opts{1} = 'NN';opts{2}='NMFinit';opts{3}='eucl';
            index = 0;
            for w_size=[0,9,24,49]
                % Looping over the window size.
                % Note that the order of chosen pixels in the window
                % is not important since there is no particular spatial
                % correlation.
                d = ones(1,R); % one element to be picked from each dictionary
                Di = cell(1,R); % generating empty multiple dictionary cell
                for r=1:R
                    Di{r} = ... % generating the dictionaries by choosing following pixels
                    Data(:,pp_index(r):min(size_synth,pp_index(r)+w_size));
                end
                % Computing the multi-dictionary model
                tic
                [~,~,K] = Dmult_nmf_OMP(Data,1000,B_init,Data(:,K_init),Di,d,opts);
                index = index+1;
                store(7+index,count,n) = length(setdiff(K(2,:),1));
                %fprintf('M2PALS-%g %d \n',w_size+1, length(setdiff(K(2,:),1)))  
                store_t(7+index,count,n) = toc;
            end
            
            % H2NMF + M2PALS
            % Automatic Dm selection by clustering R groups of pixels
            % using H2NMF.
            % This may provide poor results if multiple pure spectra are
            % lumped into one cluster. To help with this, maybe add more
            % clusters or use larger d as an upper bound.
                % Choosing the number of clusters that will serve as
                % sub-dictionaries
                Rm = R;
                % Creating the vector storing number of picked atom per
                % dictionary
                d = 2*ones(Rm,1)';
                % Computing H2NMF
                [idK, ~, ~] = hierclust2nmf(Data,Rm); 
                % Variable initialization
                Dm = cell(1,Rm);
                corr_table = cell(1,Rm);
                % Recovering Dm and the indexes matching table
                for k=1:Rm
                    % Computing Dm
                    indexes = 1:size_synth;
                    index_r = indexes(idK==k);
                    Dm{k} = Data(:,index_r);
                    % Correspondance table [1 2 3 ...; index_1 index_2 ... ]
                    corr_table{k} = [1:size(index_r,2);index_r]; 
                end
                % Computing M2PALS
                tic
                [~,~,K] = Dmult_nmf_OMP(Data,1000,B_init,Data(:,K_init),Dm,d,opts);
                % Obtaining the set of index K in the original set of
                % indexes (indexes for the DATA, not for each dictionary).
                K_true = zeros(1,R);
                for k=1:R
                    K_true(k)  = corr_table{K(1,k)}(2,K(2,k));
                end
                % Storing results
                store(12,count,n) = R - length(intersect(K_true,pp_index));
                %fprintf('H2NMF+M2PALS %d \n', R - length(intersect(K_true,pp_index)))
                store_t(12,count,n) = toc;
                
    end

end

% % Summing the matching errors over the trials mode
% store = sum(store,3);
% store_t = sum(store_t,3);
% legend = ...
% {'1-SPA','2-MPALS','3-SNPA','4-FGNSR','5-GLUP','6-NFINDR','7-SDOMP','8-M2PALS-1','9-M2PALS-10','10-M2PALS-25','11-M2PALS-50','12-H2NMF-M2PALS'};
% save('results_rev.mat','store','store_t')

% TODO : opts1 for not nonnegative case, opts4 for supress output, both with
% native default values.
% Unused material
%       FGNSR-100 for large image
%         w        =  zeros(max(IDX),1);
%         for ki = 1 : max(IDX)
%             w(ki) = sqrt(sum(IDX==ki));
%         end
        % [~,K] = fgnsr(DD*diag(w),R,ones(size(DD,2),1),'maxiter',1000,...
        % K     = K0(K);  
