function [posterior,out] = VBA_wrapup(posterior,options,dim,suffStat,u,y,it,fromPause)
% wraps up the VB inversion
% function [posterior,out] = VBA_wrapup(posterior,options,dim,suffStat,u,y,it)
% IN:
%   - [VBA internal variables; see VBA_NLStateSpaceModel.m]
% OUT:
%   - posterior: the final 'posterior' structure
%   - out: the final 'out' structure

try; fromPause; catch; fromPause = 0; end

if options.DisplayWin
    display = options.display;
end

%------- Summary output structure -----%
% Restore posterior from ODE limit if required
[posterior,options,dim,suffStat] = VBA_odeLim2NLSS(posterior,options,dim,suffStat,u,1);
% Store options, sufficient statistics, etc...
out.F = suffStat.F(end);
out.options = options;
out.u = VBA_getU(u,options,dim,'back2micro');
out.y = y;
out.dim = dim;
out.it = it;
out.suffStat = suffStat;
out.date = clock;
out.dt = toc(options.tStart);
if fromPause
    return
end


% Display inversion status
if ~options.OnLine
    [st,i] = dbstack;
    try,ifInit=isequal(st(i+2).name,'VBA_Initialize');catch,ifInit=0;end
    if ~ifInit % main VB inversion has converged
        try
            if options.DisplayWin
                % display diagnostics
                figure(display.hfp)
                VBA_ReDisplay(posterior,out);
            end
            out.options = rmfield(out.options,'display');
        end
        status = 'inversion';
    else
        status = 'initialization';
    end
    % display overall inversion time
    dt = toc(options.tStart);
    if floor(dt./60) == 0
        timeString = [num2str(floor(dt)),' sec'];
    else
        timeString = [num2str(floor(dt./60)),' min'];
    end
    str = ['VB ',status,' complete (took ~',timeString,').'];
    VBA_disp(str,options)
    VBA_disp(' ',options)
    drawnow
end




