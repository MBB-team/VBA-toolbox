function log = VBA_unit_tests()
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
log = struct('demo',{},'status',{},'stack',{},'time',{});

for i = 1:numel(demos)
   demo_name = demos{i};
   
   try
       close all
       
       fprintf('\n ####################\n  %s \n ####################\n',demo_name)
       tic ;
       evalin('base',demo_name)
       log(i) = struct('demo',demo_name,'status',1,'stack',[],'time',toc);
   catch err
       log(i) = struct('demo',demo_name,'status',0,'stack',err,'time',toc);
   end

   
   save(['../VBA_unit_test_log_' version('-release')],'log')
       
end

