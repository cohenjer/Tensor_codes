function [label, m, energy] = spkmeans(X, init)
% Perform spherical k-means clustering.
%   X: d x n data matrix
%   init: k (1 x 1) or label (1 x n, 1<=label(i)<=k) or center (d x k)
% Reference: Clustering on the Unit Hypersphere using Von Mises-Fisher Distributions.
% by A. Banerjee, I. Dhillon, J. Ghosh and S. Sra.
% Written by Michael Chen (sth4nth@gmail.com).
% Downloaded @ 
% http://www.mathworks.com/matlabcentral/fileexchange/28902-spherical-k-means/content/spkmeans.m 
% (and slightly modifed to run on previous verions of Matlab)
%% initialization
[d,n] = size(X);

if n <= init
    label = 1:init; 
    m = X; 
    energy = 0;
else

% Normalize the columns of X
X = X.*repmat( (sum(X.^2)+1e-16).^(-0.5),d,1); 

if length(init) == 1
    idx = randsample(n,init);
    m = X(:,idx);
    [ul,label] = max(m'*X,[],1);
elseif size(init,1) == 1 && size(init,2) == n
    label = init;
elseif size(init,1) == d
    m = init.*repmat( (sum(init.^2)+1e-16).^(-0.5),d,1); 
    [ul,label] = max(m'*X,[],1);
else
    error('ERROR: init is not valid.');
end
%% main algorithm: final version 
last = 0;
while any(label ~= last)
    [u,pipi,label] = unique(label);   % remove empty clusters
    k = length(u);
    E = sparse(1:n,label,1,n,k,n);
    m = X*E; 
    m = m.*repmat( (sum(m.^2)+1e-16).^(-0.5),d,1); 
    last = label;
    [val,label] = max(m'*X,[],1);
end
[ul,ul,label] = unique(label);   % remove empty clusters
energy = sum(val);
end