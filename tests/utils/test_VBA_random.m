function tests = test_VBA_random
% Unit Tests for VBA_test_random

tests = functiontests (localfunctions);

% -------------------------------------------------------------------------   
function test_arbitrary (testCase)
    testCase.verifyFail('missing test')

% -------------------------------------------------------------------------   
function test_bernoulli (testCase)

    p = rand ();

    % should fail on invalid probability
    shouldFail = @() VBA_random ('Bernoulli', -3);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');

    % should return one sample by default
    actual = VBA_random ('Bernoulli', p);
    testCase.verifyNumElements (actual, 1);
    
    % should return matrix on scalar N
    N = 3;
    actual = VBA_random ('Bernoulli', p, N);
    testCase.verifySize (actual, [N, N]);
     
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Bernoulli', p,  N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    actual = VBA_random ('Bernoulli', p, 1, 1e6);   
    % + support
    testCase.verifyEqual (sort (unique (actual)), [0 1]);
    % + mean
    testCase.verifyEqual (mean (actual), p, 'AbsTol', 1e-2);
 
% -------------------------------------------------------------------------   
function test_binomial (testCase)
    n = 2;
    p = rand ();

    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Binomial',0, p);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Binomial',n, -3);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');

    % should return one sample by default
    actual = VBA_random ('Binomial', n, p);
    testCase.verifyNumElements (actual, 1);
    
    % should return matrix on scalar N
    N = 3;
    actual = VBA_random ('Binomial', n, p, N);
    testCase.verifySize (actual, [N, N]);
     
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Binomial', n, p,  N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    actual = VBA_random ('Binomial', n, p, 1, 1e6);   
    % + support
    testCase.verifyEqual (sort (unique (actual)), 0 : n);
    % + mean
    testCase.verifyEqual (mean (actual), n * p, 'AbsTol', 1e-2);
    % + variance
    testCase.verifyEqual (var (actual), n * p * (1 - p), 'AbsTol', 1e-2);
    % + density
    for k = 0 : n
        expected = nchoosek (n, k) * (p ^ k) * ((1 - p) ^ (n - k));
        testCase.verifyEqual (mean (actual == k), expected, 'AbsTol', 1e-2);
    end  
    
% -------------------------------------------------------------------------   
function test_categorical (testCase)
    K = 3;
    p = rand (K, 1);
    p = p / sum(p);

    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Categorical', rand);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    
    % should return one sample by default
    actual = VBA_random ('Categorical', p);
    testCase.verifyNumElements (actual, 1);
    
    % should return matrix on scalar N
    N = 3;
    actual = VBA_random ('Categorical', p, N);
    testCase.verifySize (actual, [N, N]);
     
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Categorical', p,  N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    actual = VBA_random ('Categorical', p, 1, 1e6);   
    % + support
    testCase.verifyEqual (sort (unique (actual)), 1 : K);
    % + mean
    testCase.verifyEqual (mean (actual), (1 : K) * p, 'AbsTol', 1e-2);
   % + density
    for k = 1 : K
        testCase.verifyEqual (mean (actual == k), p(k), 'AbsTol', 1e-2);
    end

    % -------------------------------------------------------------------------   
function test_dirichlet (testCase)
    K = 3;
    alpha = 1 + randi (10, K, 1);

    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Dirichlet', rand);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Dirichlet', alpha, 2, 3);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    
    % should return one sample by default
    actual = VBA_random ('Dirichlet', alpha);
    testCase.verifySize (actual, [K, 1]);
    
    % should return N samples
    N = 8;
    actual = VBA_random ('Dirichlet', alpha, N);
    testCase.verifySize (actual, [K, N]);
    
    % should return sample according to distribution
    N = 1e6;
    actual = VBA_random ('Dirichlet', alpha, N);   
    % + support
    testCase.verifyTrue (VBA_isInRange (actual, [0 1]));
    testCase.verifyEqual (sum (actual), ones (1, N), 'AbsTol', 1e-13);
    % + moments
    [expected.m, expected.v] = VBA_dirichlet_moments(alpha);
    testCase.verifyEqual (mean (actual, 2), expected.m , 'AbsTol', 1e-2);
    testCase.verifyEqual (var (actual, [], 2), expected.v , 'AbsTol', 1e-2);


% -------------------------------------------------------------------------   
function test_gamma (testCase)
    a = 10 * rand ();
    b = 10 * rand ();

    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Gamma', - 1, b);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Gamma', a, - 1);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    
    % should return one sample by default
    actual = VBA_random ('Gamma', a, b);
    testCase.verifyNumElements (actual, 1);
     
    % should return matrix on scalar N
    N = 3;
    actual = VBA_random ('Gamma', a, b, N);
    testCase.verifySize (actual, [N, N]);
     
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Gamma', a, b,  N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    actual = VBA_random ('Gamma', a, b, 1, 1e6);   
    % + support
    testCase.verifyTrue (all (actual > 0));
    % + mean
    testCase.verifyEqual (mean (actual), a * b, 'RelTol', 1e-2);
    % + variance
    testCase.verifyEqual (var (actual), a * b ^ 2, 'RelTol', 1e-2);
 
  % -------------------------------------------------------------------------   
function test_gaussian (testCase)
    testCase.verifyFail('missing test')
    
% -------------------------------------------------------------------------   
function test_multinomial (testCase)
    testCase.verifyFail('missing test')
