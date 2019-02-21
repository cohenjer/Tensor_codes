function [ M ] = col_norm( M )

[m,n] = size(M);

norms = sqrt(sum(M.^2));
for r=1:n
    if norms(r)>0
        M(:,r) = M(:,r)/norms(r);
    end
end

end

