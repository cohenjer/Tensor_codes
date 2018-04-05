function tests = test_fgnsr

tests = functiontests(localfunctions);

end

function setup(testCase)

% DEBUG switch
% testCase.TestData.verbose = false;
testCase.TestData.verbose = false;

end


function test_identity_matrix(testCase)

n = 5;
M = eye(5);
r = n;

X = fgnsr(M, r, 'verbose', testCase.TestData.verbose);

testCase.assertEqual(diag(X), ones(n,1), 'AbsTol', 1e-2);

end

function test_middlepoints(testCase)

W = [
    0.0173    0.1352    0.3084
    0.3137    0.3602    0.1912
    0.3284    0.1478    0.3196
    0.1706    0.0445    0.0763
    0.1700    0.3123    0.1044
    ];

H = [
    0.5000    0.5000         0         0    1.0000         0
    0.5000         0    0.5000         0         0    1.0000
         0    0.5000    0.5000    1.0000         0         0
];

Ktrue = [4,5,6];

[~, K] = fgnsr(W*H, 3,     'verbose', testCase.TestData.verbose);
testCase.assertEqual(sort(K), sort(Ktrue));

N = randn(size(W*H));
N = 0.1 * N / norm(N, 'fro');
[~, K] = fgnsr(W*H + N, 3, 'verbose', testCase.TestData.verbose);
testCase.assertEqual(sort(K), sort(Ktrue));
    


end