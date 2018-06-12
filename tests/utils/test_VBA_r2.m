function tests = test_VBA_r2
% Unit Tests for VBA_test_r2

tests = functiontests (localfunctions);

% inputs
function test_inputs (testCase)
    shouldFail = @() VBA_r2 (rand);
    testCase.verifyError(shouldFail, '');
    shouldFail = @() VBA_r2 (rand, rand, rand, rand);
    testCase.verifyError(shouldFail, 'MATLAB:TooManyInputs');
   
% catch cases
function test_empty (testCase)
    import matlab.unittest.constraints.*
    actual = VBA_r2([], []);
    testCase.verifyThat(actual, IsScalar);
    testCase.verifyThat(actual, HasNaN);
    
% output
function test_isInRange (testCase)
    import matlab.unittest.constraints.*
    for i = 1 : 10
        pred = randn(1,1e3);
        data = 0.5 * pred;
        actual = VBA_r2(pred, data);
        testCase.verifyThat(actual, IsGreaterThanOrEqualTo(0));
        testCase.verifyThat(actual, IsLessThanOrEqualTo(1));
    end

% correct values
function test_canonical (testCase)
    N = 5e5;
    data = randn(N,1);
    for i = 5 : 10 : 100
        frac = i / 100;
        pred = data; 
        pred(1 : round (frac * N)) = mean(data);
        actual = VBA_r2(pred, data);
        testCase.verifyEqual(actual, 1-frac, 'AbsTol', 5e-3);
    end