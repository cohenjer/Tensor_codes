% From cluster indicator vector to indicator matrix

function V = vectoind(IDX,r)

m = length(IDX); 
if nargin == 1
    r = max(IDX(:)); 
end

V = zeros(m,r); 
for i = 1 : r
    V(find(IDX==i),i) = 1;    
end