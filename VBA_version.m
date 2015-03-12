function infos = VBA_version()
% VBA_VERSION
% infos = VBA_version()
% return infos about the toolbox:
%    - version: if git maintained, a branch/SHA1 tag
%    - path : local path to the toolbox

str = which('VBA_version');
infos.path = str(1:end-13);

try
    gitInf = getGitInfo ;
    infos.version = [gitInf.branch '/' gitInf.hash];
    infos.git = true;
catch
    infos.version = 'unkown (not under git control)';
    infos.git = false;
end

