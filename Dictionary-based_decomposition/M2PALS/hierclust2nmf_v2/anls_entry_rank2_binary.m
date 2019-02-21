function solve_either = anls_entry_rank2_binary(left, right)

n = size(right, 1);

solve_either = zeros(n, 2);
solve_either(:, 1) = max(0, right(:, 1) ./ left(1,1)); 
solve_either(:, 2) = max(0, right(:, 2) ./ left(2,2)); 

cosine_either = solve_either.* repmat([sqrt(left(1,1)), sqrt(left(2,2))],n,1); 

choose_first = (cosine_either(:, 1) >= cosine_either(:, 2));
solve_either(choose_first, 2) = 0;
solve_either(~choose_first, 1) = 0;