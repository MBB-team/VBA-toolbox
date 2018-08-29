function VBA_disp (str, options, carriageReturn)
% // VBA toolbox //////////////////////////////////////////////////////////
%
% VBA_disp (str, options)
% display string or cell-arry of strings according to verbose option
%
% IN:
%   - str: string or cell array of strings to display
%   - option structure containing the verbose flag (default: true)
%
% /////////////////////////////////////////////////////////////////////////

% get flag
try
    verbose = options.verbose;
catch
    verbose = true;
end

if nargin < 3
    carriageReturn = true;
end

% conditional display function
if verbose
    if iscell (str)
        for i = 1 : numel(str)
            fprintf(1, str{i});
        end
    else
       fprintf (1, str);
    end
    if carriageReturn
        fprintf ('\n');
    end
end