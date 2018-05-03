function logs = VBA_unit_tests(logs)
% This function simply launches sequentially all the demos of the toolbox and
% reports execution time and potential failures


%% find demos
vba_info = VBA_version();
[~,list]=system(['find ' vba_info.path ' -name demo_*']) ;

demos = {};
for p = strsplit(list)  
    if ~isempty(p{1})
        [~,demos{end+1},~] = fileparts(p{1});
    end
end

% setup for base
setup = 'pause off; warning off; clear all; close all; ' ;


%% run demos
logs = {}; 

parfor i = 1:numel(demos)
   demo_name = demos{i};
   
   try
       close all
       
       fprintf('\n ####################\n  %s \n ####################\n',demo_name)
       tic ;
       evalin('base',[ setup demo_name])
       logs{i} = struct('demo',demo_name,'status',1,'stack',[],'time',toc);
   catch err
       logs{i} = struct('demo',demo_name,'status',0,'stack',err,'time',toc);
   end
   
end

logs = [logs{:}];


