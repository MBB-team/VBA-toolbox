function option = check_option(option, varargin)
    if isempty(option)
        option = struct();
    end
    
    for iVar = 1:numel(varargin)/2
        varName = varargin{2*(iVar-1)+1};
        varValue = varargin{2*(iVar-1)+2};
        
        if ~isfield(option,varName)
            option.(varName) = varValue;
        end
    end
end
        
    
