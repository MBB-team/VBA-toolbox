% plot MoG approximations to the N(0,1) density
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

nmax = 4;
hf = figure('color',[1 1 1]);
na = ceil(sqrt(nmax));

for i=1:nmax
    ha(i) = subplot(na,na,i,'parent',hf);
    [m,s,w] = getMoG4N01(i,1,1);
    plot_MoG(m,s,w,struct('ha',ha(i),'gri',-5:1e-2:5));
    title(ha(i),[num2str(i),' components'])
    set(ha(i),'ylim',[0,0.5])
    grid(ha(i),'on')
end
title(ha(1),'original N(0,1) density')