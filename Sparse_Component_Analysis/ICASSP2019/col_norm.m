function [ M,N ] = col_norm( M,N )

[m,n] = size(M);

    if nargin==1
norms = sqrt(sum(M.^2));
for r=1:n
    if norms(r)>0
        M(:,r) = M(:,r)/norms(r);
    end
end
N=[];
    else
norms = sqrt(sum(M.^2));
for r=1:n
    if norms(r)>0
        M(:,r) = M(:,r)/norms(r);
        N(:,r) = N(:,r)*norms(r);
    end
end

    end
end


