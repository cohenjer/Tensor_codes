% "Fast" but less accurate SVD by computing the SVD of MM^T or M^TM 
% ***IF*** one of the dimensions of M is much smaller than the other. 
% Note. This is numerically less stable, but useful for large hyperspectral 
% images. 

function [u,s,v] = fastsvds(M,r); 

[m,n] = size(M); 
rationmn = 10; % Parameter, should be >= 1

if m < rationmn*n 
    MMt = M*M';
    [u,s,v] = svds(MMt,r); 
    v = M'*u; 
    v = v.*repmat( (sum(v.^2)+1e-16).^(-0.5),n,1); 
    s = sqrt(s); 
elseif n < rationmn*m
    MtM = M'*M;
    [u,s,v] = svds(MtM,r); 
    u = M*v; 
    u = u.*repmat( (sum(u.^2)+1e-16).^(-0.5),m,1); 
    s = sqrt(s); 
else
    [u,s,v] = svds(M,r); 
end