function tests = test_proj_omega
tests = functiontests(localfunctions);
end

function Z = cross_check(testCase, X, w, ub, diag_idx, msg)

if nargin < 6
    msg = 'Inconsistent projection results';
end

% Invoke all available algorithms and assert that results are close
Z1 = proj_omega(X, w, ub, diag_idx, 'thres');
Z2 = proj_omega(X, w, ub, diag_idx, 'parproj');

% The sparse matrix code assumes that the diagonal elements are present in
% the sparse matrix.  This condition may not be met when converting form
% dense to sparse.  We add some diagonal noise to the input matrix and hold
% thumbs.  So if this test fails, pay special attention to those
% perturbations.
if ~isempty(diag_idx)
    [m, n] = size(X);
    % Need to permute the rows of E
    E = sparse(diag_idx, 1:n, ones(1,n), m, n);
else
    E = sparse(eye(size(X)));
end

Z4 = proj_omega(X + 1e-12*E, w, ub, diag_idx, 'parproj');


% Arbitrary tolerance used, may need adjustment if this fails for too many
% false positives.
assertEqual(testCase, Z1, Z2, 'AbsTol', 1e-8, msg);
assertEqual(testCase, Z2, Z4, 'AbsTol', 1e-8, msg);

Z = Z1;

end

function test_no_inf_weight(testCase)
X = [1,0;0,1];
w = [1; +inf];

testCase.assertError(@() proj_omega(X, w, [], [],'thres'), 'FGNSR:input');

end

function test_no_neg_weight(testCase)
X = [1,0;0,1];
w = [1; -1];

testCase.assertError(@() proj_omega(X, w, [], [],'thres'), 'FGNSR:input');

end

function test_bug1(testCase)
% This turfed up a bug where the weights w were not permuted as needed.
X = [0, 1; 0, 1];
w = [1 ; 2];

Ztrue = [0, .5; 0, 1];
Ztest = cross_check(testCase, X, w, 1, []);

testCase.assertEqual(Ztrue, Ztest);


end

function test_bug2(testCase)
% Combination of zero weight and negative entries results in bound clipping
% for x(1) <= ub and x >=0.
X = [2; 2; -1];
w = [0; 1; 1];
Ztrue = [1; 2; 0];

Ztest = cross_check(testCase, X, w, 1, []);

testCase.assertEqual(Ztrue, Ztest);
end


function test_permutation1(testCase)
% Simple test for degenerated diag_index input
X = [1;1];
w = [1;2];
Ztrue = [.5; 1];
Ztest = cross_check(testCase, X, w, 1, 2);

testCase.assertEqual(Ztrue, Ztest);

end

function test_permutation2(testCase)

X = [1,0; 0,1];
diag_idx = [2,1];
Ztrue = 0.5 * ones(2,2);


% No permuatation, then X=Z
Ztest = cross_check(testCase, X, [], [], []);
testCase.assertEqual(Ztest, X);

% With permutation
Ztest = cross_check(testCase, X, [], [], diag_idx);
testCase.assertEqual(Ztest, Ztrue);


end

function test_rectangular1(testCase)
X = [1, 0; 0,1; -1, -1];

Ztrue = [1,0;0,1;0,0];
Ztest = cross_check(testCase, X, [], [], []);
testCase.assertEqual(Ztest, Ztrue);

ub = 0;
Ztrue = zeros(3,2);
Ztest = cross_check(testCase, X, [], ub, []);
testCase.assertEqual(Ztest, Ztrue);

diag_idx = [2,1];
Ztrue = [.5,.5;.5,.5;0,0];
Ztest = cross_check(testCase, X, [], [], diag_idx);
testCase.assertEqual(Ztest, Ztrue);

end

function test_random_data(testCase)
% WARNING this is a non-deterministic test

num_tries = 5;

for k=1:num_tries

seed = randi(10000000, 1, 1);
rng(seed);

m = 8;
n = 5;
X = randn(m,n);
w = rand(m,1);
ub = 3 * abs(randn(1,1));
diag_idx = randperm(n);

% Three weights 0
which = randi(m, 3, 1);
w(which) = 0.0;


msg = sprintf('Seed used: %d', seed);
cross_check(testCase, X, w, ub, diag_idx, msg);
end


end

function test_zero_weight(testCase)
X = [.5; 1.; 0.];
w = [1; 1; 0];
ub = 1.0;
Ztrue = [.75; .75; 0];
Ztest = cross_check(testCase, X, w, ub, []);
testCase.assertEqual(Ztest, Ztrue);
end
