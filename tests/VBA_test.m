function results = VBA_test (subfolder)
    vba_info = VBA_version();
    target = [vba_info.path filesep 'tests'];
    if nargin > 0
        target = [target filesep subfolder];
    end
    results = runtests(target, 'IncludeSubfolders', true);
end