function [mgp] = constrainedMoG(xt,Phi,ut,in)
% observation function for mixture of gaussian (MoG) densities
%
% function [mgp] = MoG(xt,Phi,ut,in)
% This functions serves when evaluating the best MoG approximation to a 1D
% arbitrary density.
% IN:
%   - xt: [useless]
%   - Phi: parameters of the MoG.
%   - ut: [useless]
%   - in: optinal input. This is used to pass the grid over which the MoG
%   density is evaluated.
% OUT:
%   - mgp: the MoG evaluated over the 1D grid.
%------------------------------------------------------------
% Copyright (C) 2012 Jean Daunizeau / License GNU GPL v2
%------------------------------------------------------------

try
    gri = in.gri;
catch
    gri = -10:1e-2:10;
end
gri = gri(:);

try
    cmog = in.cmog;
catch
    cmog = 1;
end


nmog = in.nmog;
s = sigm(Phi(in.s));
m = cumsum(exp(Phi(in.mu)));

if isequal(nmog/2,floor(nmog/2))
    % even #components
    if cmog ==1
        q = exp(Phi(in.q));%./2*sum(exp(Phi(in.q)));
    elseif cmog == 2
        q = exp(-0.5*m.^2.*exp(Phi(in.q)));
    end
    q = q./2.*sum(q);
    mgp = 0;
    nloops = nmog/2;
else
    % odd #components
    if cmog == 1
        q = exp(Phi(in.q));
        q0 = exp(Phi(in.q0));
    elseif cmog == 2
        q = exp(-0.5*m.^2.*exp(Phi(in.q)));
        q0 = 1;
    end
    q = q./(2*sum(q)+q0);
    q0 = q0./(2*sum(q)+q0);
    mgp = q0.*exp(-0.5*gri.^2./s)./sqrt(2*pi*s);
    nloops = (nmog-1)/2;
end

for i=1:nloops
    p1 = exp(-0.5*(gri-m(i)).^2./s);
    p1 = (1./sqrt(2*pi*s)).*p1;
    p2 = exp(-0.5*(gri+m(i)).^2./s);
    p2 = (1./sqrt(2*pi*s)).*p2;
    mgp = mgp + q(i).*p1 + q(i).*p2;
end


