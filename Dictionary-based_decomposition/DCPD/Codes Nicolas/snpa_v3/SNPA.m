function [J,H] = SNPA(M,r,normalize) 

% Successive Nonnegative Projection Algorithm (variant with f(.) = ||.||^2)
%
% *** Description ***
% At each step of the algorithm, the column of M maximizing ||.||_2 is 
% extracted, and M is updated with the residual of the projection of its 
% columns onto the convex hull of the columns extracted so far. 
% 
% See N. Gillis, Successive Nonnegative Projection Algorithm for Robust 
% Nonnegative Blind Source Separation, arXiv, 2013. 
%  
%
% [J,H] = SNPA(M,r,normalize) 
%
% ****** Input ******
% M = WH + N : a (normalized) noisy separable matrix, that is, W is full rank, 
%              H = [I,H']P where I is the identity matrix, H'>= 0 and its 
%              columns sum to at most one, P is a permutation matrix, and
%              N is sufficiently small. 
% r          : number of columns to be extracted. 
% normalize  : normalize=1 will scale the columns of M so that they sum to one,
%              hence matrix H will satisfy the assumption above for any
%              nonnegative separable matrix M. 
%              normalize=0 is the default value for which no scaling is
%              performed. For example, in hyperspectral imaging, this 
%              assumption is already satisfied and normalization is not
%              necessary. 
%
% ****** Output ******
% J        : index set of the extracted columns. 
% H        : optimal weights, that is, H argmin_{X >= 0} ||M-M(:,K)X||_F

[m,n] = size(M); 
maxitn = 100; 

if nargin <= 2, normalize = 0; end
if normalize == 1
    % Normalization of the columns of M so that they sum to one
    D = spdiags(((sum(M)+1e-16).^(-1))', 0, n, n); M = M*D; 
end

normM = sum(M.^2); 
nM = max(normM); 

i = 1; 
% Perform r recursion steps (unless the relative approximation error is 
% smaller than 10^-9)
while i <= r && max(normM)/nM > 1e-9 
    % Select the column of M with largest l2-norm
    [a,b] = max(normM); 
    % Norms of the columns of the input matrix M 
    if i == 1, normM1 = normM; end 
    % Check ties up to 1e-6 precision
    b = find((a-normM)/a <= 1e-6); 
    % In case of a tie, select column with largest norm of the input matrix M 
    if length(b) > 1, [c,d] = max(normM1(b)); b = b(d); end
    % Update the index set, and extracted column
    J(i) = b; U(:,i) = M(:,b); 
    
    % Update residual 
    if i == 1
        % Initialization using 10 iterations of coordinate descent
        H = nnlsHALSupdt(M,M(:,J),[],10); 
        % Fast gradient method for min_{y in Delta} ||M(:,i)-M(:,J)y||
        H = FGMfcnls(M,M(:,J),[],maxitn); 
    else
        H(:,J(i)) = 0; 
        h = zeros(1,n); h(J(i)) = 1; 
        H = [H; h]; 
        H = FGMfcnls(M,M(:,J),H,maxitn); 
    end
    R = M - M(:,J)*H; 
    
    % Update norms
    normM = sum(R.^2); 
   
    i = i + 1; 
end

end % of function FastSepNMF