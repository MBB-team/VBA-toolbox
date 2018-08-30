function tests = test_VBA_sigmoid
% Unit Tests for VBA_test_sigmoid

tests = functiontests (localfunctions);

% simple cases
function test_empty (testCase)
    actual = VBA_sigmoid([]);
    testCase.verifyEmpty(actual);
    
function test_canonical_scalar (testCase)
    for x = - 5 : 5
        actual = VBA_sigmoid(x);
        expected = 1 / (1 + exp (-x));
        testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
    end
    
% generalization to arbitrary dimension
function test_matrix (testCase)
    for mtx = {randn(1, 2), randn(2, 1), randn(2, 3), randn(2, 3, 4)} 
        actual = VBA_sigmoid (mtx{1});
        expected = arrayfun (@VBA_sigmoid, mtx{1});
        testCase.verifyEqual(actual, expected);
    end
    
% parametrization
function test_single_parameters (testCase)
    x = rand ();
    
    % slope
    slope = rand();
    actual = VBA_sigmoid (x, 'slope', slope);
    expected =  VBA_sigmoid (slope * x);
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
    
    % center
    center = rand();
    actual = VBA_sigmoid (x, 'center', center);
    expected =  VBA_sigmoid (x - center);
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);

    % scale
    scale = rand();
    actual = VBA_sigmoid (x, 'scale', scale);
    expected =  scale * VBA_sigmoid (x);
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
    
    % offset
    offset = rand();
    actual = VBA_sigmoid (x, 'offset', offset);
    expected =  VBA_sigmoid (x) + offset;
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
    
    % lapseRate
    shouldFail = @() VBA_sigmoid (x, 'lapseRate', - eps);
    testCase.verifyError(shouldFail, 'MATLAB:InputParser:ArgumentFailedValidation');
    shouldFail = @() VBA_sigmoid (x, 'lapseRate', 0.5 + eps);
    testCase.verifyError(shouldFail, 'MATLAB:InputParser:ArgumentFailedValidation');
    
    rate = 0.5 * rand ();
    actual = VBA_sigmoid (x, 'lapseRate', rate);
    expected =  VBA_sigmoid (x, 'offset', rate, 'scale', 1 - 2 * rate) ;
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
 
% inverse
function test_inverse (testCase)
    x = rand ();

    actual = VBA_sigmoid (VBA_sigmoid (x, 'inverse', true));
    expected =  x;
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-9);
    
% derivatives
function test_derivatives (testCase)
    x = rand (1, 2);
    
    % no arg
    [~, actual.dsdx, actual.dsdp] = VBA_sigmoid (x);
    expected =  diag(VBA_numericDiff(@VBA_sigmoid, 1, x))';
    testCase.verifyEqual(actual.dsdx, expected, 'AbsTol', 1e-5);
    testCase.verifyEmpty(actual.dsdp);
    
    % parametrized
    parameters = {'slope', 'center', 'scale', 'offset'}; 
    for iP = 1 :numel (parameters) 
        p = 0.5 * rand();
        [~, ~, actual] = VBA_sigmoid (x, parameters{iP}, p);
        expected =  VBA_numericDiff(@VBA_sigmoid, 3, x, parameters{iP}, p);
        testCase.verifyEqual(actual, expected, 'AbsTol', 1e-5);
    end
 
    % multiple parameters
    p = rand (1, 2);
    [~, ~, actual] = VBA_sigmoid (x, 'slope', p(1), 'center', p(2));
    [~, ~, expected(1,:)] = VBA_sigmoid (x, 'slope', p(1), 'center', p(2), 'derivatives',{'center'});
    [~, ~, expected(2,:)] = VBA_sigmoid (x, 'slope', p(1), 'center', p(2), 'derivatives',{'slope'});
    testCase.verifyEqual(actual, expected, 'AbsTol', 1e-5);

    
