function ha = plot_MoG(m,s,w,in)

% plots 1D mixture of Gaussian densities
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

try; type = get(in.ha,'type'); catch; type = []; end
if isequal(type,'axes')
    ha = in.ha;
else
    hf = figure('color',[1 1 1]);
    ha = axes('parent',hf);
end
set(ha,'nextplot','add');
xlabel('x')
ylabel('MoG weighted components p(x)')

try
    gri = in.gri;
catch
    gri = -20:1e-2:20;
end
gri = gri(:);

nmog = length(w);
mgp = 0;
for i=1:nmog
    ppi = exp(-0.5*(gri-m(i)).^2./s(i))./sqrt(2*pi*s(i));
    mgp = mgp + w(i).*ppi;
    plot(ha,gri,w(i).*ppi,'r--')
end
plot(ha,gri,mgp,'k');
