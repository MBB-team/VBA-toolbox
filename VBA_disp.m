function VBA_disp(str,options)

% conditional display function
if options.verbose
    if iscell(str)
        n = length(str);
        for i=1:n
            disp(str{i})
        end
    else
       disp(str) 
    end
end