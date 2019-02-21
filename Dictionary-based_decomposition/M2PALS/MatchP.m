function [indices] = MatchP(X,D)

%[~,m] = size(X);
%[~,K] = size(D);

Gram     =    X'*D;
%norms    =    sqrt(sum(X.^2));

[~,indices] = max(Gram'); 
    
%    for k=1:K
%        
%        inner_ik  =     Gram(i,k);%/norms(i);%/norm(D(:,i)); 
%        if inner_ik>inner
%           inner    =    inner_ik;
%           indice   =    k;
%        end
%        
%    end
%    


   
end