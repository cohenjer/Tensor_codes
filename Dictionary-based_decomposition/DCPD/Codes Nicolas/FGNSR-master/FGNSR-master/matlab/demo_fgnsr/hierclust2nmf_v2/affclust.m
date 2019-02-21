% Display clusters

function [a, Vaff] = affclust(K,H,L,ncol,bw); 

n = H*L;

if iscell(K)
    K = clu2vec(K,n); 
end

% K is an indicator vector of type IDX
A = vectoind(K); 
r = size(A,2); 

% 'Optimize' display in 16/9
if nargin < 4
    ncol = 1; nrow = ceil(r/ncol);
    while (r > 1 && L*ncol*9 < H*nrow*16) || rem(r,ncol) == 1
        ncol = ncol+1; 
        nrow = ceil(r/ncol);
    end
end
if nargin == 5 && bw == 1
    [a, Vaff] = affichage(A,ncol,H,L,bw); %bw==1 switches black and white
else
    [a, Vaff] = affichage(A,ncol,H,L); 
end