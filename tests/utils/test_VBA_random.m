function tests = test_VBA_random
% Unit Tests for VBA_random

tests = functiontests (localfunctions);

% -------------------------------------------------------------------------   
function test_arbitrary (testCase)
    K = 3;
    p = rand (K, 1);
    p = p / sum (p);
    vals_1 = 'abc';
    vals_k = ['ab'; 'cd'; 'ef'];
    M = size(vals_k, 2);
    
    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Arbitrary', [0.5; 0.5], vals_1);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Arbitrary', zeros (3, 1), vals_1);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Arbitrary', p, vals_k, 8, 8);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    
    % should return one sample by default
    actual = VBA_random ('Arbitrary', p, vals_1);
    testCase.verifyNumElements (actual, 1);
    actual = VBA_random ('Arbitrary', p, vals_k);
    testCase.verifySize (actual, [M, 1]);
    
    % should return matrix on scalar N
    N = 8;
    actual = VBA_random ('Arbitrary', p, vals_1, N);
    testCase.verifySize (actual, [N, N]);
    actual = VBA_random ('Arbitrary', p, vals_k, N);
    testCase.verifySize (actual, [M, N]);
    
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Arbitrary', p, vals_1, N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    % + univariate
    actual = VBA_random ('Arbitrary', p, vals_1, 1, 1e6);   
    % - support
    testCase.verifyEqual (unique (actual), vals_1);
    % + density
    for i = 1 : numel (p)
        testCase.verifyEqual (mean (actual == vals_1(i)), p(i), 'AbsTol', 1e-2);
    end 
    % + multivariate
    actual = VBA_random ('Arbitrary', p, vals_k, 1e6);   
    % - support
    testCase.verifyEqual (unique (actual', 'rows'), vals_k);
    % + density
    for i = 1 : numel (p)
        isEq = all(bsxfun(@eq, actual, vals_k(i,:)'));
        testCase.verifyEqual (mean (isEq, 2), p(i), 'AbsTol', 1e-2);
    end 

% -------------------------------------------------------------------------   
function test_bernoulli (testCase)
    p = rand ();

    % should fail on invalid probability
    shouldFail = @() VBA_random ('Bernoulli', -3);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    % should fail on invalid size
    shouldFail = @() VBA_random ('Bernoulli', rand (3), 2);
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
    
    % should return matrix on matrix p
    actual = VBA_random ('Bernoulli', rand (N{:}));
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    actual = VBA_random ('Bernoulli', p, 1, 1e6);   
    % + support
    testCase.verifyEqual (unique (actual), [0 1]);
    % + mean
    testCase.verifyEqual (mean (actual), p, 'AbsTol', 1e-2);
    
    % should deal with nan
    actual = VBA_random ('Bernoulli', [p NaN]);   
    testCase.verifyTrue (~ isnan (actual(1)));
    testCase.verifyTrue (isnan (actual(2)));
    
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
    testCase.verifyEqual (unique (actual), 0 : n);

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
    shouldFail = @() VBA_random ('Categorical', zeros (K, 1));
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
    testCase.verifyEqual (unique (actual), 1 : K);
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
    
    % + univariate
    mu_1 = randn ();
    Sigma_1 = rand () ^ 2;
    % + multivariate
    k = 2;
    mu_k = randn (k, 1);
    Sigma_k = rand (k);
    Sigma_k = Sigma_k * Sigma_k';
    
     % should fail on invalid parameters
    shouldFail = @() VBA_random ('Gaussian', mu_1, Sigma_k);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Gaussian', mu_k, Sigma_k, 2, 2);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    
    % should return one sample by default
    % + univariate
    actual = VBA_random ('Gaussian', mu_1, Sigma_1);
    testCase.verifyNumElements (actual, 1);
    % + multivariate
    actual = VBA_random ('Gaussian', mu_k, Sigma_k);
    testCase.verifySize (actual, [k, 1]);

    % should return matrix on scalar N
    N = 5;
    % + univariate
    actual = VBA_random ('Gaussian', mu_1, Sigma_1, N);
    testCase.verifySize (actual, [N, N]);
    % + multivariate
    actual = VBA_random ('Gaussian', mu_k, Sigma_k, N);
    testCase.verifySize (actual, [k, N]);
     
    % should return matrix of asked dimension
    N = num2cell (1 + randi (5, 1, 4));
    actual = VBA_random ('Gaussian', mu_1, Sigma_1,  N{:});
    testCase.verifySize (actual, [N{:}]);
    
    % should return sample according to distribution
    % + univariate
    actual = VBA_random ('Gaussian', mu_1, Sigma_1, 1, 1e6);   
    % - mean
    testCase.verifyEqual (mean (actual), mu_1, 'AbsTol', 1e-2);
    % - variance
    testCase.verifyEqual (var (actual), Sigma_1, 'AbsTol', 1e-2);
    % + multivariate
    actual = VBA_random ('Gaussian', mu_k, Sigma_k, 1e6);   
    % - mean
    testCase.verifyEqual (mean (actual, 2), mu_k, 'AbsTol', 1e-2);
    % - variance
    testCase.verifyEqual (cov (actual'), Sigma_k, 'AbsTol', 1e-2);
  
% -------------------------------------------------------------------------   
function test_multinomial (testCase)
    n = 2;
    K = 3;
    p = rand (K, 1);
    p = p / sum(p);
    
    % should fail on invalid parameters
    shouldFail = @() VBA_random ('Multinomial',0, p);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Multinomial',n, .5);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Multinomial',n, [1; 1]);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');
    shouldFail = @() VBA_random ('Multinomial',n, p, 8, 8);
    testCase.verifyError(shouldFail, 'VBA:invalidInput');

    % should return one sample by default
    actual = VBA_random ('Multinomial', n, p);
    testCase.verifySize (actual, [K, 1]);
    
    % should return matrix on scalar N
    N = 8;
    actual = VBA_random ('Multinomial', n, p, N);
    testCase.verifySize (actual, [K, N]);
    
    % should return sample according to distribution
    N = 1e6;
    actual = VBA_random ('Multinomial', n, p, N);   
    % + support
    testCase.verifyEqual (unique (actual)', 0 : n);
    testCase.verifyEqual (sum (actual), n * ones (1, N), 'AbsTol', 1e-13);
    % + density
    for i = 1 : K
        expected = n * p(i);
        testCase.verifyEqual (mean (actual(i, :)), expected, 'AbsTol', 1e-2);
        expected = n * p(i) * (1 - p(i));
        testCase.verifyEqual (var (actual(i, :)), expected, 'AbsTol', 1e-2);
    end  
    
    % should deal with nan
    actual = VBA_random ('Multinomial', n, nan(K, 1));
    testCase.verifySize (actual, [K, 1]);
    testCase.verifyTrue (all (isnan (actual)));
    