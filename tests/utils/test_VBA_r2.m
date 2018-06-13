function tests = test_VBA_r2
% Unit Tests for VBA_test_r2

tests = functiontests (localfunctions);

% inputs
function test_inputs (testCase)
    shouldFail = @() VBA_r2 (rand);
    testCase.verifyError(shouldFail, 'MATLAB:narginchk:notEnoughInputs');
    shouldFail = @() VBA_r2 (rand, rand, rand, rand);
    testCase.verifyError(shouldFail, 'MATLAB:TooManyInputs');
   
% catch overly simple cases
function test_empty (testCase)
    import matlab.unittest.constraints.*
    for k = 0 : 2
        actual = VBA_r2 (ones (1, k), ones (1, k));
        testCase.verifyEqual (actual, nan);
    end
    
% output
function test_isInRange (testCase)
    data = randn (1, 1e2);
    actual = VBA_r2 (data, data);
    testCase.verifyEqual (actual, 1);
    actual = VBA_r2 (- data, data);
    testCase.verifyEqual (actual, 0);

% correct values
function test_canonical (testCase)
    N = 5e5;
    data = randn (N, 1);
    for i = 5 : 20 : 100
        frac = i / 100;
        pred = data; 
        pred(1 : round (frac * N)) = mean (data);
        actual = VBA_r2 (pred, data);
        testCase.verifyEqual (actual, 1-frac, 'AbsTol', 5e-3);
    end
    
% data exclusion
function test_exclusion (testCase)
    actual = VBA_r2 ([0 1 0], [0 1 1], [0 0 1]);
    testCase.verifyEqual (actual, 1);
