
% display evaluation of VB inversion (missing region pb)
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

lo.theta = [1];
lo.sigma = [1e0 1e1 1e2 1e4];
TR = [.5 1 2 4];
N = 20;

display = 'stochastic'; % 'deterministic', or 'stochastic'


for MR = 1:2
    hff(MR) = figure('color',ones(1,3));
    if MR ==1
        str = ' : 1 --> (2) --> 3';
    else
        str = ' : (2) --> 1 and 3';
    end
    set(hff(MR),'name',['missing region: pb #',num2str(MR),str])
    for i=1:length(lo.sigma)
        for j=1:length(TR)
            ha = subplot(4,4,(i-1)*4+j,'parent',hff(MR));
            dAC = mr(MR).dAC{i,j};
            dA = reshape(dAC(1:4),2,2);
            sAC = mr(MR).sAC{i,j};
            sA = reshape(sAC(1:4),2,2);
            
            if isequal(display,'deterministic')
                imagesc(dA,'parent',ha)
            else
                imagesc(sA,'parent',ha)
            end
            set(ha,'clim',[-1 1])
            title(ha,...
                ['SNR=',num2str(lo.sigma(i)),...
                ' ; TR=',num2str(TR(j))])
            box(ha,'on')
            grid(ha,'on')
        end
    end
end

