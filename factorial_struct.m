function structList = factorial_struct(varargin)
% FACTORIAL_STRUCT 
% structList = factorial_struct('param1', {value1a, value1b,...}, 'param2', 1:n, ...)
% provide a structure with fields 'param1', 'param2', ... whose values are
% set in a factorial design.
% eg:
%
% structList(1).param1 = value1_a;
% structList(1).param2 = 1;
%
% structList(2).param1 = value1_a;
% structList(2).param2 = 2;
% ...
%
% structList(n).param1 = value1_a;
% structList(n).param2 = n;
%
% structList(n+1).param1 = value1_b;
% structList(n+1).param2 = 1;
% ...
%
% structList(2*n).param1 = value1_b;
% structList(2*n).param2 = n;
%
    
    % ________________________________________________________
    % initialize struct
    structList = struct();

    % ________________________________________________________
    % recusrive field adjonction
    if ~isempty(varargin)
 
        % extract first variable name and potential values
    	argName      = varargin{1};
        argValueList = varargin{2};
        
        if ~iscell(argValueList)
            try, argValueList = num2cell(argValueList); end
        end
        
        % get structure for remaining parameters
        subStructList = factorial_struct(varargin{3:end}); 
       
        % iterate over potential value, 
        for iValue = 1:numel(argValueList)
            [subStructList.(argName)] = deal(argValueList{iValue});
            structList = struct_cat(structList,subStructList);
        end
        
        % sort field by calling order
        structList = orderfields(structList,varargin(1:2:end));
    end
    
end

function s = struct_cat(s1,s2)
    try
        s = [s1 , s2];
    catch
        s=s2;
    end
end
