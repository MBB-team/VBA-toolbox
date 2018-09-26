function tests = test_VBA_exceedanceProbability
% Unit Tests for VBA_exceedanceProbability

tests = functiontests (localfunctions);

% -------------------------------------------------------------------------   
function test_empty (testCase)
    actual = VBA_exceedanceProbability ([]);
    testCase.verifyEmpty (actual);
    actual = VBA_exceedanceProbability ([], []);
    testCase.verifyEmpty (actual);