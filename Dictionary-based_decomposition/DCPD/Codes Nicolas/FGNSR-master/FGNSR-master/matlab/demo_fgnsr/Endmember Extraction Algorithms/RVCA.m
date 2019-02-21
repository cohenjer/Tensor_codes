% Robust VCA 

function K = RVCA(M,r,rparam); 

maxiter = 10; 
if nargin <= 2
    rparam = 10;
end

emin = +Inf; 

for i = 1 : rparam
    [A, K] = VCA(M,'Endmembers',r,'verbose','off');  
    H = nnlsHALSupdt(M,M(:,K),[],maxiter);
    err = norm(M-M(:,K)*H,'fro'); 
    
    if  err < emin
        Kf = K; 
        emin = err;
    end
end