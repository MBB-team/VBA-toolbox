function results = VBA_test (subfolder)
    import matlab.unittest.TestSuite;        

    vba_info = VBA_version();
    target = [vba_info.path filesep 'tests'];
    if nargin > 0
        target = [target filesep subfolder];
    end
    suite = TestSuite.fromFolder(target, 'IncludingSubfolders', true);
    results = suite.run ();
end