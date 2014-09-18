function VBA_disp(str,options)

% conditional display function
if options.verbose
    if iscell(str)
        n = length(str);
        for i=1:n
            fprintf(str{i})
        end
    else
       fprintf(str) 
    end
    fprintf('\n')
end