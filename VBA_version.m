function [infos, status] = VBA_version()
% infos = VBA_version()
% return infos about the current version of the toolbox. It requires the
% toolbox to be installed using git (not by downlading the zip) to work properly:
%       git clone git@github.com:MBB-team/VBA-toolbox.git
% OUT: 
%   - infos.path: path of the toolbox in use
%   * under git versioning, the function also returns:
%   - infos.version: the current branch and revision of the toolbox (branch/commit) 
%   - infos.status: a summary of the status of the current revision compared to
%              the one found on github.

str = which('VBA_version');
infos.path = str(1:end-13);

try % if the toolbox was installed through git
    gitInf = getGitInfo ;
    infos.version = [gitInf.branch '/' gitInf.hash];
    infos.git = true;
    
    if nargout == 2
        
    try % try to see if new commits are online
        request = sprintf('https://api.github.com/repos/MBB-team/VBA-toolbox/compare/%s...%s',gitInf.branch,gitInf.hash);
        tracker=webread(request);
        
        switch tracker.status
            case 'identical'
                status = 'The toolbox is up to date.';
            case 'behind'
                status = sprintf('The toolbox is %d revision(s) behind the online version.',tracker.behind_by);
            case 'ahead'
                status = sprintf('The toolbox is %d revision(s) ahead the online version.',tracker.ahead_by);
        end

    end
    end
        
catch
    infos.version = 'unkown (not under git control)';
    infos.git = false;
end

