function V = make_perm(R,k,n)
% This function gives all sets of indexes from n:R of cardinal k, without
% repeting the permuted index sets.
V=[];

if k>R
    V = [];
elseif k==1
    V = linspace(n,R,R-n+1)';
elseif k==R
    V = linspace(n,R,R-n+1);
else
    for i=n:R-k+1
    V = [V;i*ones(nchoosek(R-i,k-1),1),make_perm(R,k-1,i+1)];    
    end
end

end