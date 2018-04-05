function example_midpoints

disp('Running a small example...');
disp('==========================');


% The vertices of the convex hull are at positions 4,6,7,8, all other
% points are middle points between these vertices

H = [
    0.5  0    0.5  0    0    1.0  0    0    0    0.5
    0.5  0.5  0    0    0    0    1.0  0    0.5  0
    0    0.5  0.5  0    0.5  0    0    1.0  0    0
    0    0    0    1.0  0.5  0    0    0    0.5  0.5
];


% Generate conic basis at random
m = 6;
n = size(H, 2);
r = 4;
W = rand(m,r);
D = diag(1./sum(W));
W = W*D;

% Compute data matrix
M = W*H;

% Generate noise of norm 0.01
N = rand(m,n);
N = 0.01 * N / norm(N, 'fro');

% Recover the basic indices using FGNSR
[X, K] = fgnsr(M + N, r, 'maxiter', 100);
disp('True basic indices: 4,6,7,8');
disp('Computed indices:')
disp(sort(K(:))');
disp('Diagonal of X:')
disp(diag(X)');

disp('==========================');
disp('done.');
end