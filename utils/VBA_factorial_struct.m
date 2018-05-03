function structList = VBA_factorial_struct (varargin)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% structList = VBA_factorial_struct (param1Name, param1Values, param2Name, param2Values, ...)
% Create an array of structures with fields 'param1Name', 'param2Name', ... 
% whose values, across across array elements, will cover the cartesian
% product of param1Values x param2Values x ...
%
% IN:
%   - list of pairs of arguments:
%       - a string representing the label
%       - an array or cell-array defining the potential values for the
%         label
% OUT:
%   - array of structure defining the factorial space of input values, ie,
%     all possible combinations of label/values
%
% Example
% ~~~~~~~
%
% >> design = VBA_factorial_struct( ...
%   'condition', {'A', 'B'}, ... 
%   'intensity', 0 : 10 : 10 );
%
% will create an 4 x 1 array with the following structures
% design(1):
%    condition: 'A'
%    intensity: 0
% design(2):
%    condition: 'A'
%    intensity: 10
% design(3):
%    condition: 'B'
%    intensity: 0
% design(4):
%    condition: 'B'
%    intensity: 10
%
% >> struct2table(design)
%
% struct2table(design)
% ans =
%   6×2 table
%     condition    intensity
%     _________    _________
%     'A'           0       
%     'A'          10       
%     'B'           0       
%     'B'          10       
%     
% /////////////////////////////////////////////////////////////////////////
    
% initialize array
% =========================================================================
structList = [];

% endpoint of recursion
% =========================================================================
if isempty(varargin)
    return
end
 
% recursive field adjonction
% =========================================================================
% extract first variable name and potential values
argName      = varargin{1};
argValueList = varargin{2};

% values should be a cell array
if ~ iscell (argValueList)
    try
        argValueList = num2cell (argValueList); 
    end
end
        
% get structure for remaining parameters recursively
subStructList = VBA_factorial_struct(varargin{3:end}); 
       
% iterate over values and add the field
for iValue = 1 : numel (argValueList)
    [subStructList.(argName)] = deal (argValueList{iValue});
    structList = [structList ; subStructList];
end
        
% sort field by calling order
structList = orderfields (structList, varargin(1:2:end));
    
end

