function [] = VBA_pause(options)
% used to pause the VB inversion from the GUI interactively
% NB: when paused, the inversion allows interactive diagnosis...

if ~ options.DisplayWin
    return;
end
try
    dt = toc(options.tStart);
    timeString = sprintf('Elapsed time: %d min and %d sec', floor(dt/60), round(rem(dt,60)));
    set(options.display.htt, 'string', timeString);
end

try
    hpause = options.display.hpause;
    if ~isempty(hpause) && ishandle(hpause)
        if get(hpause,'value')
            set(hpause,'string','PAUSED!',...
                'backgroundColor',[1 0.5 0.5])
            stop = 0;
            try
                [posterior,out] = evalin(...
                    'caller',...
                    'VBA_wrapup(posterior,options,options.dim,suffStat,suffStat.u,y,it,1)');
            catch
                [posterior,out] = evalin(...
                    'caller',...
                    'VBA_wrapup(posterior,options,options.dim,suffStat,u,y,it,1)');
            end
            hfp = VBA_ReDisplay(posterior,out,1,1);
            s = dbstatus;
            if isempty(s)
                dbstop if error % this allows to have a look...
            end
            while ~stop
                pause(2)
                if ~ ishandle(hfp)
                    set(hpause, 'value', 0);
                end
                if ~ get(hpause,'value') 
                    stop = 1;
                    set(hpause,'string','pause and diagnose?',...
                        'backgroundColor',0.8*[1 1 1])
                    try
                        close(hfp)
                    end
                end
            end
        end
    end
end


