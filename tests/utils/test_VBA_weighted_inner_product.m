function tests = test_VBA_weighted_inner_product
% Unit Tests for VBA_weighted_inner_product

tests = functiontests (localfunctions);

function test_input_checks (testCase)
% Test input parameter checking

shouldFail = @() VBA_weighted_inner_product(zeros(2,3), zeros(3,3));
testCase.verifyError(shouldFail, 'VBA:invalidInput');
shouldFail = @() VBA_weighted_inner_product(zeros(3,3), zeros(3,3), zeros(4,1));
testCase.verifyError(shouldFail, 'VBA:invalidInput');
shouldFail = @() VBA_weighted_inner_product(zeros(4,3), zeros(4,3), zeros(2,2));
testCase.verifyError(shouldFail, 'VBA:invalidInput');
shouldFail = @() VBA_weighted_inner_product(zeros(3,3), zeros(3,3), zeros(4,5));
testCase.verifyError(shouldFail, 'VBA:invalidInput');

function test_weight_vector (testCase)
% Test basic functioning

m = 10;
n_A = 3;
n_B = 4;
range = [0, 1000]; % test with random integers, so equality testing can be exact
A = randi(range, m, n_A);
B = randi(range, m, n_B);
w = randi(range, m, 1);

actual = VBA_weighted_inner_product(A, B, w);
testCase.verifySize(actual, [n_A, n_B]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), w.*B(:,3)))

actual = VBA_weighted_inner_product(A, B, w');
testCase.verifySize(actual, [n_A, n_B]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), w.*B(:,3)))


function test_weight_matrix (testCase)
% Test basic functioning

m = 10;
n_A = 3;
n_B = 4;
range = [0, 1000]; % test with random integers, so equality testing can be exact
A = randi(range, m, n_A);
B = randi(range, m, n_B);
W = randi(range, m, m);

actual = VBA_weighted_inner_product(A, B, W);
testCase.verifySize(actual, [n_A, n_B]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), W*B(:,3)))

function test_optional_args (testCase)
% Test optional arguments

m = 10;
n_A = 3;
n_B = 4;
range = [0, 1000]; % test with random integers, so equality testing can be exact
A = randi(range, m, n_A);
B = randi(range, m, n_B);
w = randi(range, m, 1);

actual = VBA_weighted_inner_product(A, [], w);
testCase.verifySize(actual, [n_A, n_A]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), w.*A(:,3)))

actual = VBA_weighted_inner_product(A, [], []);
testCase.verifySize(actual, [n_A, n_A]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), A(:,3)))

actual = VBA_weighted_inner_product(A);
testCase.verifySize(actual, [n_A, n_A]);
testCase.verifyEqual(actual(2,3), dot(A(:,2), A(:,3)))
