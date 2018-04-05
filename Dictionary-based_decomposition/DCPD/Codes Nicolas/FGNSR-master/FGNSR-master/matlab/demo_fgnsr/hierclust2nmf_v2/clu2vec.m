% Transform a cluster cell to a vector 

function IDX = clu2vec(K,m,r); 

if nargin < 3
    r = length(K);
end
if nargin < 2
    % Compute max entry in K
    m = 0; 
    for i = 1 : r
        m = max(0, max(K{i})); 
    end
end
IDX = zeros(m,1); 
for i = 1 : r 
    IDX(K{i}) = i; 
end