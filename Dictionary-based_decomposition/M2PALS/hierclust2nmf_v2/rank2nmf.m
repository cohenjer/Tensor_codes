% Given a data matrix M (m-by-n), computes a rank-two NMF of M. 
%
% See Algorithm 3 in 
% 
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv.  
% 
% ****** Input ******
%  M     : a nonnegative m-by-n data matrix  
%
% ****** Output ******
%  (U,V) : a rank-two NMF of M
%  s1    : first singular value of M 

function [U,V,s1] = rank2nmf(M)

[m,n] = size(M); 

% Best rank-two approximation of M
if min(m,n) == 1
    [U,S,V] = fastsvds(M,1); 
    U = abs(U); V = abs(V); s1 = S; 
else
    [u,s,v] = fastsvds(M,2); 
    s1 = s(1); 
    K = FastSepNMF(s*v',2); 
    U = zeros(size(M,1),2); 
    if length(K) >= 1
        U(:,1) = max(u*s*v(K(1),:)',0); 
    end
    if length(K) >= 2
        U(:,2) = max(u*s*v(K(2),:)',0); 
    end
    % Compute corresponding optimal V 
    V = anls_entry_rank2_precompute_opt(U'*U, M'*U); 
end