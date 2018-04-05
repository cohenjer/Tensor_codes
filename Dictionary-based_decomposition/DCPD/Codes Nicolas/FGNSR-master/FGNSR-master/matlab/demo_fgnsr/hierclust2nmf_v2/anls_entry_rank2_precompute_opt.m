% Solve min_H ||M - WH'||_2 s.t. H >= 0
%
% where left = W^TW and right = M^TW 
% 
% See Kuang, Park, `Fast Rank-2 Nonnegative Matrix Factorization 
% for Hierarchical Document Clustering', KDD '13. 
%
% See also Algorithm 4 in 
% Gillis, Kuang, Park, `Hierarchical Clustering of Hyperspectral Images 
% using Rank-Two Nonnegative Matrix Factorization', arXiv. 
%
% ****** Input ******
%  left     : 2-by-2 matrix (or possibly 1-by-1)
%  right    : n-by-2 matrix (or possibly n-by-1)
%
% ****** Output ******
%   H    : nonnegative n-by-2 matrix, solution to KKT equations

function H = anls_entry_rank2_precompute_opt(left, right)

warning('off'); 
if length(left) == 1
    H = max(0,right/left);
else
    H = (left \ right')'; 
    use_either = ~all(H>=0, 2);
    H(use_either, :) = anls_entry_rank2_binary(left, right(use_either,:)); 
    H = H'; 
end
warning('on'); 


% Case where one entry in each column of H has to be equal to zero
function solve_either = anls_entry_rank2_binary(left, right)

n = size(right, 1);

solve_either = zeros(n, 2);
solve_either(:, 1) = max(0, right(:, 1) ./ left(1,1)); 
solve_either(:, 2) = max(0, right(:, 2) ./ left(2,2)); 

cosine_either = solve_either.* repmat([sqrt(left(1,1)), sqrt(left(2,2))],n,1); 

choose_first = (cosine_either(:, 1) >= cosine_either(:, 2));
solve_either(choose_first, 2) = 0;
solve_either(~choose_first, 1) = 0;