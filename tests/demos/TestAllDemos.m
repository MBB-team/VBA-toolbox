classdef TestAllDemos < matlab.unittest.TestCase

    properties (TestParameter)
      demo = findDemos ();
    end
       
    methods (Test)
        function testDemo (testCase, demo)         
            setup = 'pause off; warning off; ' ;
            try
                evalin ('base', [setup demo ' ; close all'])
                %evalc ([setup demo ' ; close all']);
            catch err
                testCase.verifyFail (err.message);
            end
        end
    end
end

function demos = findDemos ()
    
    vba_info = VBA_version();
    d = dir(fullfile(vba_info.path, '**', 'demo_*.m'));
    
    demos = {};
    for i = 1 : numel (d)  
        demos{end+1} = d(i).name(1:end-2);
    end
    
end