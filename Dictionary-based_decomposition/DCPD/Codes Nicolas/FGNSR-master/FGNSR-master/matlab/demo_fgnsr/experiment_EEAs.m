% Running the algorithm on data set M 

function [Kall, Hall, resultsErr, resultsTime, J, nEEAs, IDX] = experiment_EEAs( M , r , nbsubs )

% 1. Identify the subsampled data matrix using H2NMF 
fprintf('Identifying the %3.0f centroid pixels... \n',nbsubs); 
[IDX,~,J] = hierclust2nmf(M,nbsubs); 
MS = M(:,J); % Subsampled data matrix

% 2. Running the EEAs on the full and subsampled data set
maxitNNLS = 500; 
nEEAs{1} = 'SPA  '; 
nEEAs{2} = 'VCA  '; 
nEEAs{3} = 'XRAY '; 
nEEAs{4} = 'H2NMF'; 
nEEAs{5} = 'SNPA '; 
for algo = 1 : 5
    for s = 1 : 2
        % Dictionary D 
        if s == 1
            D = M;   % full data set
            fprintf([nEEAs{algo}, ' on the full data set...']); 
        else
            D = MS; % subsampled data set
            fprintf(' on the centroids...');
        end
        %% Normalization to apply EEAs that need it (SPA, VCA, H2NMF, SNPA) 
        %D = D./repmat(sum(D),size(D,1),1); 
        % EEA
        e = cputime;
        K = EEAs(D,r,algo); % extract endmembers from 'dictionary'
        resultsTime(algo+(s-1)*5) = cputime - e;
        Kall{algo+(s-1)*5} = K;
        % Error
        H = nnlsHALSupdt(M,D(:,K),[],maxitNNLS); % optimal weights
        Hall{algo+(s-1)*5} = H;
        resultsErr(algo+(s-1)*5) = 100*norm(M - D(:,K)*H,'fro')/norm(M,'fro');
    end
    fprintf(' done.\n'); 
end

% 3. Running NSRFG on the subsampled data set 
fprintf('FGNSR on the centroids...'); 
algo = 11;    
maxiterFG = 1000;  
% Weigthing the columns
for ki = 1 : max(IDX)
    w(ki) = sqrt(sum(IDX==ki));
end

[X, K] = fgnsr(MS*diag(w), r, ones(size(MS,2),1), 'maxiter', maxiterFG, 'verbose', false, 'proj', 'parproj', 'col_weights', ones(size(MS,2),1) ); 

resultsTime(algo) = cputime - e;
Kall{algo} = K;
% Error
H = nnlsHALSupdt(M,MS(:,K),[],maxitNNLS);
Hall{algo} = H;
resultsErr(algo) = 100*norm(M - MS(:,K)*H,'fro')/norm(M,'fro');
disp('done.')