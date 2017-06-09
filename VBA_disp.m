function VBA_disp(str,options)
try
    verbose = options.verbose;
catch
    verbose = 1;
end
% conditional display function
if verbose
    if iscell(str)
        n = length(str);
        for i=1:n
            fprintf(1,str{i})
        end
    else
       fprintf(1,str) 
    end
    fprintf('\n')
end