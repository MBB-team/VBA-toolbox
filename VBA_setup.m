function VBA_setup()
% VBA_setup()
% Cleanly add the VBA_toolbox to the Matlab path.
% Please prefer this function to addpath() or addpath(genpath()) to ensure
% all subfolders are included with the notable exception of the .git tracker.

fprintf('\n');
fprintf(' ======================================================================================\n')
fprintf(' VBA-toolbox installation \n')
fprintf(' ======================================================================================\n')
fprintf('\n');


%%

isWin = ~isempty(strfind(computer, 'PCWIN')) || strcmp(computer, 'i686-pc-mingw32');
isOSX = ~isempty(strfind(computer, 'MAC')) || ~isempty(strfind(computer, 'apple-darwin'));
isLinux = strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64') || ~isempty(strfind(computer, 'linux-gnu'));


%%
% Locate ourselves:
targetdirectory=fileparts(mfilename('fullpath'));
if ~strcmpi(targetdirectory, pwd)
    fprintf(' *** You need to change your working directory to the VBA-toolbox folder before running this routine!\n');
    return
else
    fprintf(' Will setup %s as the working copy of the VBA-toolbox folder.\n',targetdirectory);
    
end
fprintf('\n');


%%
% Does SAVEPATH work?
err=savepath;

if err
    % warn user the new path wont be saved:
    p=fullfile(matlabroot,'toolbox','local','pathdef.m');
    
    fprintf(' *** savepath() failed. Probably the pathdef.m file lacks write permission. \n\n');
    
    fprintf('     You will have to re-run %s again after each Matlab startup untill you fix this permission issue.\n',mfilename);
    
    fprintf('     Please ask a user with administrator privileges to enable write by everyone for the file:\n')
    fprintf('     ''%s''\n\n',p);
    
    fprintf(' // Installation will continue but your changes will not be saved.\n')
    
    
    
end


%%
% Remove old "VBA-toolbox" from path:
if isOSX || isLinux
    p = strsplit(path,':');
elseif isWin
    p = strsplit(path,';');
end
    
r = regexp(p, '^.*VBA-toolbox.*$','match');
r = r(~cell2mat(cellfun(@isempty,r,'UniformOutput',false)));

if ~isempty(r)
    try
        % try to get version
        ver = VBA_version();
    catch
        % find root path
        pl = cell2mat(cellfun(@(x) numel(x{1}),r,'UniformOutput',false));
        ridx = find(pl == min(pl));
        
        ver.path = r{ridx};
        ver.version = 'unknown';
    end
    
    
    fprintf(' *** The VBA-toolbox already appears in your path:\n');
    fprintf('     - installation path: %s\n',ver.path);
    fprintf('     - release: %s\n',ver.version);
    
    
    fprintf('     Installation will remove this toolbox from your path (without deleting the files).\n');
    
    cont = NaN;
    while isnan(cont)
        answer = input('\n     Are you sure you want to continue (yes/no)? ','s');
        switch answer
            case 'yes'
                fprintf('\n //  Alright, let''s go...\n')
                cont = true;
            case 'no'
                fprintf('\n //  You''re the boss! \n');
                cont = false;
            otherwise
                fprintf('\n //  I didn''t understand, please try again.\n');
        end
    end
    
    if ~cont
        fprintf('\n Run %s again if you change your mind.\n\n',mfilename);
        return
    end
    
    % remove all old path
    
    for iP = 1:numel(r)
        rmpath(str2mat(r{iP}))
    end
    try
        savepath
    end
end

%% find folder to install
p = strsplit(genpath(pwd),':');
ref_path = ['refs' filesep 'heads'];
r = regexp(p, ['^.*VBA-toolbox' filesep '.git' filesep '(?!' ref_path ').*$'],'match');
p = p(cell2mat(cellfun(@isempty,r,'UniformOutput',false)));

for iP = 1:numel(p)
    addpath(str2mat(p{iP}))
end
try
    savepath
end

fprintf('\n')
try
    ver = VBA_version();
    fprintf(' The VBA-toolbox has been successfully installed in %s\n\n',ver.path);
    fprintf(' Type demo_Qlearning to give it a try\n\n');
    fprintf(' ======================================================================================\n');
    
    
catch err
    fprintf(' *** Something wrong happened:\n')
    fprintf('     -  %s \n\n',err.message);
end



end


