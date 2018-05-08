function [x,gx,microTime,sampleInd] = VBA_microTime(posterior,u,out)
% gets micro-time time series from posterior
try
    inF = out.options.inF{1};
catch
    inF = out.options.inF;
end
options = out.options;
try
    dt = options.decim.*inF.deltat;
    microTime = 0:inF.deltat:options.dim.n_t*dt;
catch
    microTime = 0:options.dim.n_t*options.decim;
end
sampleInd = 1+[1:options.dim.n_t].*options.decim;

% pre-allocate variables
x = zeros(options.dim.n,options.dim.n_t*options.decim+1);
gx = zeros(options.dim.p,options.dim.n_t*options.decim+1);

% deal with empty input
if isempty(u)
    u = zeros(0,options.dim.n_t);
end

% Initial conditions
x(:,1) = posterior.muX0;
gx(:,1) = NaN; % no observation on initial condition!

%--- Loop over time series ---%
for t=1:options.dim.n_t
    if options.dim.n_theta >= 1
        if options.OnLine
            theta = posterior.muTheta(:,t);
        else
            theta = posterior.muTheta;
        end
    else
        theta = [];
    end
    if options.dim.n_phi >= 1
        if options.OnLine
            phi = posterior.muPhi(:,t);
        else
            phi= posterior.muPhi;
        end
    else
        phi = [];
    end
    for i=1:options.decim
        indT = (t-1).*options.decim+i+1;
        if options.microU
            tu = indT-1;
        else
            tu = t;
        end
        if options.skipf(t)
            x(:,indT) = f_Id(x(:,indT-1),theta,u(:,tu),[]);
        else
            x(:,indT) = options.f_fname(x(:,indT-1),theta,u(:,tu),options.inF);
        end
        if isfield(out,'suffStat') && isfield(out.suffStat,'dx') && ...
                ~isempty(out.suffStat.dx) && isequal(i,options.decim)
            % add state noise to evolution mapping
            x(:,indT) = x(:,indT) + out.suffStat.dx(:,t);
        end
        [gx(:,indT)] = options.g_fname(x(:,indT),phi,u(:,tu),options.inG);
    end
end



