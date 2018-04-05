% Extract "most" representative column from a matrix M as follows: 
% 
% First, it computes the best rank-one approximation u v^T of M. 
% Then, it identifies the column of M minimizing the MRSA with the first
% singular vector u of M. 
% 
% See Section 4.4.1 of 
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 

function [u,s,b] = reprvec(M); 

[u,s,v] = svds(M,1); 
u = abs(u); 
[m,n] = size(M); 
% Exctract the column of M approximating u the best (up to a translation and scaling)
u = u - mean(u); 
Mm = M - repmat(mean(M),m,1); 
err = acos(  (Mm'*u/norm(u))./( sqrt(sum(Mm.^2)) )' ); 
[a,b] = min(err); 
u = M(:,b); 