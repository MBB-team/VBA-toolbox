function results = VBA_test ()
    vba_info = VBA_version();
    results.utils = runtests([vba_info.path filesep 'tests'], 'IncludeSubfolders', true, 'UseParallel', true);
end