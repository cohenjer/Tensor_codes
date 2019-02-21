% Different EEA algorithms

function K = EEAs(M,r,algo); 

if algo == 1 
    K = FastSepNMF(M,r); 
elseif algo == 2 
    K = VCA(M,'Endmembers',r,'verbose','off');  
elseif algo == 3
    K = FastConicalHull(M,r); 
elseif algo == 4
    [~, ~, K] = hierclust2nmf(M,r,1,[],0);
elseif algo == 5
    K = SNPA(M,r); 
end