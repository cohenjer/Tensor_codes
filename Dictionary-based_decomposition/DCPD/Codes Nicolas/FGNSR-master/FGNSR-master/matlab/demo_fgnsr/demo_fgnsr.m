function demo_fgnsr

% Compare endmember extraction algorithms (EEAs) on the subsampled Urban
% hyperspecral image (every second pixel in each direction kept).
% 
% The EEAs under consideration are 
% VCA, SPA, XRAY, SNPA, H2NMF, and FGNSR. 
% 
% See also: Gillis & Luce, "A Fast Gradient Method for Nonnegative Sparse
% Regression with Self Dictionary"

% 0. Set up FGNSR and other path
old_pwd = pwd;
cd('fgnsr-0.1');
setup_fgnsr();
cd(old_pwd);

addpath('hierclust2nmf_v2'); 
addpath('Endmember Extraction Algorithms');  

% 1. Load (undersampled) Urban data set into 'M'
load('urban_undersampled4.mat', 'M');

% Image size
dimim = 77;
% Number of pure pixels
r = 6; 
% Number of subsampled pixels to be able to run NSRFG 
nbsubs = 100; 

% 2. Run experiment 
[Kall, Hall, resultsErr, resultsTime, J, nEEAs] = experiment_EEAs( M , r , nbsubs ) ; 

% 3. Display results 
fprintf('------------------------------------------------------ \n')
fprintf('  Algorithm   | Total time (s.) |   Relative Error (%%) \n')
fprintf('------------------------------------------------------ \n')
for i = 1 : 5
    fprintf([nEEAs{i},'         |       %3.2f      |        %2.2f  \n'], resultsTime(i), resultsErr(i));
end
for i = 1 : 5
    fprintf([nEEAs{i}, '-%3.0f' , '     |       %3.2f      |        %2.2f  \n'], nbsubs, resultsTime(i+5), resultsErr(i+5));
end
fprintf('FGNSR-%3.0f     |       %3.2f      |        %2.2f  \n', nbsubs, resultsTime(11), resultsErr(11));
fprintf('------------------------------------------------------ \n')
% Display results of NSRFG
figure; 
plot( M(:, J(Kall{11})), '-'); 
title('Spectral signatures extrated by FGNSR'); 
[~,n] = size(M); 
show_abmap(Hall{i}',r,dimim,n/dimim); 
title('Abundance maps identified by FGNSR'); 


end
