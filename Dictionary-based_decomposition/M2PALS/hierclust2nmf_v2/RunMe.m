% Example of using H2NMF on the Urban dataset

close all; clear all; clc; 

load Urban; 
% Select a subset of wavelengths (17) to reduce computational time
% (Remove this line to have the same test as in the paper)
R = R(:,:,1:10:end); 

% Run H2NMF
algo = 1; r = []; 
[IDX, C, J, sol] = hierclust2nmf(R,r,algo); 

% Display clusters
close all; nimrow = 3; 
affclust(IDX,307,307,nimrow);  
title('Clusters','FontSize',18); 

% Display endmembers
figure; plot(C); 
title('Endmembers','FontSize',18);  
axis([1 size(C,1) 0 max(C(:))*1.05]); 

% Compute and display abundance maps of endmembers
for i = 1 : size(R,3);
    M(i,:) = reshape(R(:,:,i),307*307,1);
end 
H = nnlsm_blockpivot(C,M); 
affichage(H'.^(0.75),nimrow,307,307); 
title('Abundance Maps','FontSize',18); 