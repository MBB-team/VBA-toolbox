function option = VBA_check_struct(option, varargin)
% compare and fill-in input structure with default structure

if isempty(option)
    option = struct();
end

if numel(varargin) == 1 && isempty(varargin{1})
    % if nothing given, do nothin
    return;
elseif numel(varargin) == 1 && isstruct(varargin{1})
    % default struct as template
    varName = fieldnames(varargin{1});
    varValue = struct2cell(varargin{1});
else
    % default as a list of key values pair
    assert(mod(numel(varargin),2)==0,'*** default values should be defined as a structure or a cell of key value pairs.') ;
    for iVar = 1:numel(varargin)/2
        varName{iVar} = varargin{2*(iVar-1)+1};
        varValue{iVar} = varargin{2*(iVar-1)+2};
    end
end

for iVar = 1:numel(varName)
    if ~isfield(option,varName{iVar})
        option.(varName{iVar}) = varValue{iVar};
    end
end


