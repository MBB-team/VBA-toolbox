function VBA_disp(str,options)

% conditional display function

if options.verbose
    disp(str)
end