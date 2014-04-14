function [X0,phi,theta] = VBA_multisessionUnique(options)

if isfield(options,'multisession') ...
        && isfield(options.multisession,'expanded') ...
        && options.multisession.expanded
    multisession = options.multisession;
    
    uniX0 = setdiff(multisession.indices.X0(:,1),multisession.fixed.X0);
    X0 = multisession.indices.X0(uniX0,2:end);
    X0 = [multisession.indices.X0(:,1) ; X0(:)] ;
    
    uniPhi = setdiff(multisession.indices.phi(:,1),multisession.fixed.phi);
    phi = multisession.indices.phi(uniPhi,2:end);
    phi = [multisession.indices.phi(:,1) ; phi(:)] ;
    
    uniTheta = setdiff(multisession.indices.theta(:,1),multisession.fixed.theta);
    theta = multisession.indices.theta(uniTheta,2:end);
    theta = [multisession.indices.theta(:,1) ; theta(:)] ;
    
    if isequal(options.g_fname,@VBA_odeLim)
        dim = options.inG.old.dim;
        theta = theta + dim.n_phi;
        X0 = X0 + (dim.n_phi + dim.n_theta);
        phi = [phi ; theta; X0];
        theta = [];
        X0 = [];
        
    end
else
    X0 = 1:options.dim.n;
    phi = 1:options.dim.n_phi;
    theta = 1:options.dim.n_theta;
end

% 

end