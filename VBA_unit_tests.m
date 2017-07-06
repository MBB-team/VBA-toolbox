function log = VBA_unit_tests()
% This function simply launches sequentially all the demos of the toolbox and
% reports execution time and potential failures

d = dir('./subfunctions/demo*');

log = struct('demo',{},'status',{},'stack',{},'time',{});

for i = 1:numel(d)
   demo_name = d(i).name(1:end-2);
   
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

