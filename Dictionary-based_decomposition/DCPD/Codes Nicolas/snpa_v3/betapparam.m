% Given an m-by-r matrix W, this function comptues the parameter beta'(W)
%
% See Remark 1 in the SNPA paper for more details. 
%
% Reordering the columns of W in the order they are extracted by SNPA, 
% 
% beta'(W) = betap = the minimum over 1 <= i < r between 
%      min_{j > i} ||R_{W(:,1:i)} W(:,j)||_2 and 
%      min_{i < j < k} ||R_{W(:,1:i)} W(:,j) - R_{W(:,1:i)} W(:,k)||_2 ) ] 
% 
% where  R_{W(:,1:i)} W(:,j) = W(:,j) - W(:,1:i)x, 
%            with x = argmin_{x>=0, sum_i x_i = 1} ||W(:,j)-W(:,1:i)x||_2
%
% Note that it is rather slow to compute (e.g., ~3s. for a 20 x 20 matrix). 

function betap = betapparam(W) 

[m,r] = size(W); 
% Order columns of W as extracted by SNPA 
K = SNPA(W,r);  
W = W(:,K); 
betap =  minnormdist(W); 
for i = 1 : r-1
    % compute residuals R_{W(:,i)}(W(:,i+1:r))
    Yi = FGMfcnls(W(:,i+1:r),W(:,1:i)); 
    Wi = W(:,i+1:r) - W(:,1:i)*Yi; 
    betap =  min(betap,minnormdist(Wi)); 
end


% Given a matrix W, comptues 
%
% minnd = min( min_j ||W(:,j)||_2 , min_{i ~= j} ||W(:,i) - W(:,j)||_2 )

function minnd = minnormdist(W) 

[m,r] = size(W); 
minnd =  min(sqrt(sum(W.^2))); 
for i = 1 : r
    for j = i+1 : r
        minnd = min( minnd , norm(W(:,i) - W(:,j))/sqrt(2) ); 
    end
end