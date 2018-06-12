function results = VBA_test ()
    vba_info = VBA_version();
    results = runtests([vba_info.path filesep 'tests'], 'IncludeSubfolders', true,'UseParallel', true);
end