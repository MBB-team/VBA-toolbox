function tests = test_VBA_exceedanceProbability
% Unit Tests for VBA_exceedanceProbability

tests = functiontests (localfunctions);

% general checks
% =========================================================================
function test_fails_on_unknown_distribution (testCase)
    actual = @() VBA_exceedanceProbability ('ugassnai', [1 1], eye(2));
    testCase.verifyError(actual, '')
    
% gaussian case
% =========================================================================

% input cheks
% -------------------------------------------------------------------------   
function test_gaussian_fails_on_empty_moment (testCase)
    actual = @() VBA_exceedanceProbability ('Gaussian', [], []);
    testCase.verifyError (actual, '')
    
function test_gaussian_return_nan_for_univariate (testCase)
    actual = VBA_exceedanceProbability ('Gaussian', 1, 1);
    testCase.verifyEqual (actual, NaN)
    
% computations
% -------------------------------------------------------------------------     
function test_gaussian_unity_vector (testCase)
    % ep should be the same size as mu and sum to unity
    k = 5;
    mu = rand(1,k);
    Sigma = diag(randn(1,k).^2);
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
    testCase.verifySize (actual, [k, 1]);
    testCase.verifyEqual (sum (actual), 1, 'AbsTol', eps);

function test_gaussian_bivariate (testCase)
    % bivariate case, known analöytical solution
    k = 2;
    mu = randn(1, k);
    Sigma = diag(randn(1,k).^2);
    % cdf of the difference of the two varaibles
    expected = .5 * erfc (((mu(2) - mu(1)) / sqrt (2 * sum (diag( Sigma)))));
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
    testCase.verifyEqual (actual(1), expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma, struct('method','sampling'));
    testCase.verifyEqual (actual(1), expected, 'AbsTol', 1e-3);
    
function test_gaussian_centered (testCase)
    % centered case, equal eps
    k = 5;
    mu = zeros(1, k);
    Sigma = randn()^2 * eye(k);
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3);
    
function test_gaussian_infinite_precision (testCase)
    % centered case and infinite precision, equal eps
    k = 5;
    mu = zeros(1, k);
    Sigma = zeros(k);
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3);

function test_gaussian_limit_case (testCase)
    % limit case, indicator ep
    k = 5;
    mu = [1, zeros(1, k - 1)];
    Sigma = zeros(k);
    expected = [1 ; zeros(k - 1, 1)];
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3);

 function test_gaussian_sampling_approx_analytical (testCase)
     % general case, sampling and analytical should give the same answer
     k = 5;
     mu = rand(1, k);
     Sigma = diag(randn(1, k).^2);
     actual = VBA_exceedanceProbability ('Gaussian', mu, Sigma, struct('method','sampling'));
     expected = VBA_exceedanceProbability ('Gaussian', mu, Sigma);
     testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3);
    
% dirichlet case
% =========================================================================

% input cheks
% -------------------------------------------------------------------------   
function test_dirichlet_fails_on_empty_moment (testCase)
    actual = @() VBA_exceedanceProbability ('Dirichlet', []);
    testCase.verifyError (actual, '')
    
function test_dirichlet_return_nan_for_univariate (testCase)
    actual = VBA_exceedanceProbability ('Dirichlet', 1);
    testCase.verifyEqual (actual, NaN);
    
function dirichlet_fails_on_negative (testCase)
    actual = @() VBA_exceedanceProbability ('Dirichlet', [-1 0]);
    testCase.verifyError (actual, '')
    
% computations
% -------------------------------------------------------------------------     
function test_dirichlet_unity_vector (testCase)
    % ep should be the same size as alpha and sum to unity
    k = 5;
    alpha = randi(5,1,k);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifySize (actual, [k, 1]);
    testCase.verifyEqual (sum (actual), 1, 'AbsTol', eps);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifySize (actual, [k, 1]);
    testCase.verifyEqual (sum (actual), 1, 'AbsTol', eps);
    
function test_dirichlet_bivariate (testCase)
    % bivariate case, analytical solution
    k = 2;
    alpha = randi(5, 1, k);
    expected = betainc (0.5, alpha(1) , alpha(2), 'upper');
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual(1), expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual(1), expected, 'AbsTol', 1e-3); 
    
function test_dirichlet_centered (testCase)
    % equifrequent case, equal eps
    k = 5;
    alpha = k * ones(1, k);
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3); 
  
function test_dirichlet_centered_small_alphas (testCase)
    % equifrequent case, equal eps
    k = 5;
    alpha = ones(1, k) / k;
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3); 
    
  function test_dirichlet_infinite_precision (testCase)
    % centered case and infinite precision, equal eps
    k = 5;
    alpha = 1e4 * ones(1, k);
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3);
  
  function test_dirichlet_infinite_variance (testCase)
    % centered case and infinite precision, equal eps
    k = 5;
    alpha = 1e-4 * ones(1, k);
    expected = ones (k, 1) / k;
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3); 

function test_dirichlet_limit_case (testCase)
    % limit case, indicator ep
    k = 5;
    alpha = [Inf, eps * ones(1, k - 1)];
    expected = [1 ; zeros(k - 1, 1)];
    actual = VBA_exceedanceProbability ('Dirichlet', alpha);
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3); 

 function test_dirichlet_sampling_approx_analytical (testCase)
     % general case, sampling and analytical should give the same answer
     k = 5;
     alpha = randi(5,1,k);
     actual = VBA_exceedanceProbability ('Dirichlet', alpha, struct('method','sampling'));
     expected = VBA_exceedanceProbability ('Dirichlet', alpha);
     testCase.verifyEqual (actual, expected, 'AbsTol', 1e-3); 
     