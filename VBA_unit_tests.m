function logs = VBA_unit_tests(logs)
% This function simply launches sequentially all the demos of the toolbox and
% reports execution time and potential failures


%% find demos
[~,list]=system('find . -name demo*') ;

demos = {};
for p = strsplit(list)  
    if ~isempty(p{1})
        [~,demos{end+1},~] = fileparts(p{1});
    end
end


%% run demos
logs = {}; 

parfor i = 1:numel(demos)
   demo_name = demos{i};
   
   try
       close all
       
       fprintf('\n ####################\n  %s \n ####################\n',demo_name)
       tic ;
       evalin('base',demo_name)
       logs{i} = struct('demo',demo_name,'status',1,'stack',[],'time',toc);
   catch err
       logs{i} = struct('demo',demo_name,'status',0,'stack',err,'time',toc);
   end
   
end

logs = [logs{:}];


