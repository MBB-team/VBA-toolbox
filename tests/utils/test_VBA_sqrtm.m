function tests = test_VBA_sqrtm
% Unit Tests for VBA_sqrtm

tests = functiontests (localfunctions);

% -------------------------------------------------------------------------   
function test_empty (testCase)
    actual = VBA_sqrtm ([]);
    testCase.verifyEmpty (actual);
    
function test_scalar (testCase)
    C = 4;
    actual = VBA_sqrtm (C);
    testCase.verifyEqual (actual, sqrt (C));
    actual = VBA_sqrtm (C, true);
    testCase.verifyEqual (actual, 1 / sqrt (C));
    
function test_zero (testCase)
    C = zeros(4);
    actual = VBA_sqrtm (C);
    testCase.verifyEqual (actual, C);
    actual = @() VBA_sqrtm (C, true);
    testCase.verifyError (actual, 'VBA:expectedFinite');
    
function test_infinite (testCase)
    C = diag (inf(1, 3));
    % infinite diagonal
    actual = VBA_sqrtm (C);
    testCase.verifyEqual (actual, C);
    % infinite precision
    actual = VBA_sqrtm (C, true);
    testCase.verifyEqual (actual, zeros (size (C)) );
    % ill defined
    actual = @() VBA_sqrtm (inf(3, 3));
    testCase.verifyError (actual, 'VBA:expectedFinite'); 
    
 function test_diagonal (testCase)
    C = diag (rand (1, 3));
    C = C * C';
    actual = VBA_sqrtm (C);
    testCase.verifyEqual (actual, sqrt (C));
    actual = VBA_sqrtm (C, true);
    expected = diag (1 ./ sqrt ( diag (C)));
    testCase.verifyEqual (actual, expected, 'AbsTol', 1e-11);
    
 function test_general (testCase)
     C = randn (3);
     C = C * C';
     actual = VBA_sqrtm (C);
     testCase.verifyEqual (actual * actual, C, 'AbsTol', 1e-11);
     actual = VBA_sqrtm (C, true);
     testCase.verifyEqual (inv (actual * actual), C, 'AbsTol', 1e-11);


 
  
    