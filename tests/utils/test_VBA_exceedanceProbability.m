function tests = test_VBA_exceedanceProbability
% Unit Tests for VBA_exceedanceProbability

tests = functiontests (localfunctions);

% gaussian case
% -------------------------------------------------------------------------   
function test_gaussian_empty (testCase)
    actual = VBA_exceedanceProbability ([], []);
    testCase.verifyEmpty (actual);
  
function test_gaussian_unity (testCase)
    % ep for univariate should be one
    mu = rand();
    Sigma = rand();
    actual = VBA_exceedanceProbability (mu, Sigma);
    testCase.verifyEqual (actual, 1);
    
    % ep should be the same size as mu and sum to unity
    k = 5;
    mu = rand(1,k);
    Sigma = diag(randn(1,k).^2);
    actual = VBA_exceedanceProbability (mu, Sigma);
    testCase.verifyNumElements (actual, k);
    testCase.verifyEqual (sum (actual), 1, 'AbsTol', eps);

function test_gaussian_canonical (testCase)
    % bivariate case, analytical solution
    k = 2;
    mu = randn(1, k);
    Sigma = randn(1,k).^2;
    actual = VBA_exceedanceProbability (mu, diag(Sigma));
    expected = .5 * erfc (((mu(2) - mu(1)) / sqrt (2 * sum (Sigma)))); % cdf diff
    testCase.verifyEqual (actual(1), expected, 'AbsTol', 1e-4);
    % centered case, equal eps
    k = 5;
    mu = zeros(1, k);
    Sigma = diag(randn(1,k).^2);
    actual = VBA_exceedanceProbability (mu, Sigma);
    expected = ones (k, 1) / k;
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-12);
    
 function test_gaussian_sampling (testCase)
     options.method = 'sampling';
     % empty
     actual = VBA_exceedanceProbability ([], [], options);
     testCase.verifyEmpty (actual);
     % unity
     actual = VBA_exceedanceProbability (1, 1, options);
     testCase.verifyEqual (actual, 1);
     % general case
     k = 2;
     mu = rand(1,k);
     Sigma = diag(randn(1,k).^2);
     actual = VBA_exceedanceProbability (mu, Sigma, options);
     expected = VBA_exceedanceProbability (mu, Sigma);
     testCase.verifyLessThan (norm (actual - expected), 1e-3);
    
% dirichlet case
% -------------------------------------------------------------------------   
function test_dirichlet_empty (testCase)
    actual = VBA_exceedanceProbability ([]);
    testCase.verifyEmpty (actual);
    
    
function test_dirichlet_unity (testCase)
    mu = rand();
    Sigma = rand();
    actual = VBA_exceedanceProbability (mu, Sigma);
    testCase.verifyEqual (actual, 1);