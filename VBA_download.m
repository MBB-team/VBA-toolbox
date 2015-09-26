function VBA_download()


IsWin = ~isempty(strfind(computer, 'PCWIN')) || strcmp(computer, 'i686-pc-mingw32');
IsOSX = ~isempty(strfind(computer, 'MAC')) || ~isempty(strfind(computer, 'apple-darwin'));
IsLinux = strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64') || ~isempty(strfind(computer, 'linux-gnu'));


%% check if under git control
infos = VBA_version();

fprintf('\n')

if ~infos.git
    fprintf('  This repository is not under git control, and cannot be updated automatically.\n')
    fprintf('  You need to reinstall the VBA-toolbox as following:\n')
    fprintf('   - Go to the folder where you want to install the toolbox;\n')
    fprintf('   - Clone the latest toolbox using in your terminal (<a href="http://git-scm.com/downloads">help</a> for installing git):\n\n')
    fprintf('     git clone https://github.com/MBB-team/VBA-toolbox.git\n\n')
    fprintf('   - Go the the newly created VBA-toolbox folder;\n')
    fprintf('   - Reinstall the toolbox using VBA_setup() in Matlab;\n')
    fprintf('   - Delete the old VBA-toolbox folder.\n')
    fprintf('\n')
    
    return
end

%% check if git is installed
gitpath = GetGitPath;

if isempty(gitpath)
   fprintf('\n')
   fprintf('  Could not locate git.\n')
   fprintf('  Please make sure git is correctly installed (<a href="http://git-scm.com/downloads">help</a> for installing git).\n\n')
end

%% move to toolbox directory
toolbox_dir = fileparts(mfilename('fullpath'));
origin_dir = pwd;
cd(toolbox_dir);



%% find current head
%git rev-parse  HEAD
getHeadCommand=[gitpath 'git rev-parse  HEAD'];

if IsOSX || IsLinux
    [err, result]=system(getHeadCommand);
   % result = 'For reason, see output above.';
else
    % MK: TODO: Check if this is still needed on >= R2007a on Windows?
    % I think it was only for R11 on Windows.
    [err, result]=dos(getHeadCommand, '-echo');
end

if err
    fprintf('  Sorry. The update command failed:\n');
    fprintf('  %s\n', result);
    return
end

if strcmp(['master\' result], infos.version)
    fprintf('  The toolbox is already up to date\n\n')
    return
else
    fprintf('  Preparing to dowlad revision %s \n', result)
    
   % git stash
   % git fetch
   % git merge
   % git stash pop

   % getUpdateCommand=[gitpath 'git fetch']

end



cd(origin_dir);



    