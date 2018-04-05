function H = nchoose2(r)
% H = nchoose2(r)
%
% Generate a matrix H r-by-N, with N = r(r-1)/2 such that each column of H
% has exactly two non-zero elements equal to 1 in different positions
% 
% ****** Input ******
% r : positive integer 
% ****** Output ******
% H : an r-by-N matrix as described above

R = nchoosek(r,2);
H = zeros(r,R);
k = 1; 
for i = 1 : r
    for j = i+1:r
        H(i,k) = 1; H(j,k) = 1; 
        k = k + 1; 
    end
end

end % of function nchoose2