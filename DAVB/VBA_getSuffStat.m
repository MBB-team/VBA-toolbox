function [suffStat] = VBA_getSuffStat(options,suffStat,flag)
% fills in default sufficient statistics structure
% function [suffStat] = VBA_getSuffStat(options)
% IN:
%   - options: the options structure (only the substructure .dim is used
%   here)
%   - suffStat: {[]}
%   - flag: {0}, 1 if online version

try
    flag;
catch
    flag = 0;
end

try
    suffStat.gx;
catch
    suffStat.gx = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.vy;
catch
    suffStat.vy = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.muX;
catch
    suffStat.muX = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.dy;
catch
    suffStat.dy = zeros(options.dim.p,options.dim.n_t);
end
try
    suffStat.dx;
catch
    suffStat.dx = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.vdx;
catch
    suffStat.vdx = zeros(options.dim.n,options.dim.n_t);
end
try
    suffStat.dx0;
catch
    suffStat.dx0 = zeros(options.dim.n,1);
end
try
    suffStat.dy2;
catch
    suffStat.dy2 = 0;
end
try
    suffStat.dx2;
catch
    suffStat.dx2 = 0;
end
try
    suffStat.dtheta;
catch
    if ~flag
        suffStat.dtheta = zeros(options.dim.n_theta,1);
    else
        suffStat.dtheta = zeros(options.dim.n_theta,options.dim.n_t);
    end
end
try
    suffStat.dphi;
catch
    if ~flag
        suffStat.dphi = zeros(options.dim.n_phi,1);
    else
        suffStat.dphi = zeros(options.dim.n_phi,options.dim.n_t);
    end
end
try
    suffStat.Ssigma;
catch
    suffStat.Ssigma = 0;
end
try
    suffStat.Sphid2gdphi2;
catch
    suffStat.Sphid2gdphi2 = 0;
end
try
    suffStat.Sphi;
catch
    suffStat.Sphi = 0;
end
try
    suffStat.Sphid2gdphidx;
catch
    suffStat.Sphid2gdphidx = 0;
end
try
    suffStat.Sthetad2fdtheta2;
catch
    suffStat.Sthetad2fdtheta2 = 0;
end
try
    suffStat.Stheta;
catch
    suffStat.Stheta = 0;
end
try
    suffStat.Sthetad2fdthetadx;
catch
    suffStat.Sthetad2fdthetadx = 0;
end
try
    suffStat.Salpha;
catch
    suffStat.Salpha = 0;
end
try
    suffStat.SXd2fdx2;
catch
    suffStat.SXd2fdx2 = 0;
end
try
    suffStat.SXtdfdx;
catch
    suffStat.SXtdfdx = 0;
end
try
    suffStat.SXd2gdx2;
catch
    suffStat.SXd2gdx2 = 0;
end
try
    suffStat.trSx;
catch
    suffStat.trSx = 0;
end
try
    suffStat.SX;
catch
    suffStat.SX = 0;
end
try
    suffStat.SX0;
catch
    suffStat.SX0 = 0;
end

if options.binomial
    try
        suffStat.logL;
    catch
        suffStat.logL = -Inf;
    end
end
