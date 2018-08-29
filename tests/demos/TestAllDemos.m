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
    [~,list]=system(['find ' vba_info.path filesep 'demos -name demo_*.m']) ;

    demos = {};
    for p = strsplit(list)  
        if ~isempty(p{1})
            [~,demos{end+1},~] = fileparts(p{1});
        end
    end
    
end