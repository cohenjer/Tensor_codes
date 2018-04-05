% Given an m-by-r matrix W, comptues:
%
% alpha = min_i min_{x >= 0} ||W(:,i) - W(:,I) x||_2 / ||W(:,i)||_2
% where I = {1,2,...,r}\{i}
% 
% and k = index for which the minimum is achieved. 

function [alpha,k] = alphaparam(W) 

[m,r] = size(W); 
% Normalization
for i = 1 : r 
    W(:,i) = W(:,i)/norm(W(:,i),2);
end
alpha = +Inf;
% Looping over the column of W
for i = 1 : r
    x = nnlsHALSupdt(W(:,i),W(:,[1:i-1 i+1:r]),[],100); 
    alphai =  norm(W(:,i) - W(:,[1:i-1 i+1:r])*x ); 
    if alphai < alpha
        alpha = alphai;
        k = i;
    end
end