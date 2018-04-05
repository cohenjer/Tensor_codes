% Given a matrix M, split its columns into two subsets
% 
% See Section 3 in 
%
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 
%
%
% ****** Input ******
%  M     : m-by-n data matrix (or a H-by-L-by-m tensor) 
%  algo  : algorithm used to split the clusters
%          1. rank-two NMF (default)
%          2. k-means
%          3. spherical k-means
%
% ****** Output ******
%   K    : two clusters 
%   U    : corresponding centroids
%   s    : first singular value of M

function [K,U,s] = splitclust(M,algo); 

if nargin == 1
    algo = 1;
end
if algo == 1  % rank-2 NMF
    [U,V,s] = rank2nmf(M); 
    % Normalize columns of V to sum to one
    V = V.*repmat( (sum(V)+1e-16).^(-1), 2,1); 
    x = V(1,:)'; 
    % Compute treshold to split cluster 
    threshold = fquad(x); 
    K{1} = find(x >= threshold); 
    K{2} = find(x < threshold);  
    
elseif algo == 2 % k-means
    [u,s,v] = fastsvds(M,2); % Initialization: SVD+SPA
    Kf = FastSepNMF(s*v',2,0);
    U0 = u*s*v(Kf,:)'; 

    [IDX,U] = kmeans(M', 2, 'EmptyAction','singleton','Start',U0'); 
    U = U'; 
    K{1} = find(IDX==1); 
    K{2} = find(IDX==2); 
    s = s(1); 
    
elseif algo == 3 % shperical k-means
    [u,s,v] = fastsvds(M,2); % Initialization: SVD+SPA 
    Kf = FastSepNMF(s*v',2,0);
    U0 = u*s*v(Kf,:)'; 
    
    [IDX,U] = spkmeans(M, U0); 
    % or (?)
    %[IDX,U] = kmeans(M', 2, 'EmptyAction','singleton','Start',U0','Distance','cosine'); 
    K{1} = find(IDX==1); 
    K{2} = find(IDX==2); 
    s = s(1); 
end