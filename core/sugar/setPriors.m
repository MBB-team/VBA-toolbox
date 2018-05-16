function [dim, options] = setPriors(options, dim, priorType,varargin)

    switch priorType
        case 'theta'
            [options.priors.muTheta,options.priors.SigmaTheta, options.inF.thetaLabel] = feval(@priorUglyfier, varargin{:});  
            nTheta=length(options.priors.muTheta);
            dim.n_theta = nTheta;
        case 'phi'
            [options.priors.muPhi,options.priors.SigmaPhi, options.inG.phiLabel] =  feval(@priorUglyfier, varargin{:});
            nPhi = length(options.priors.muPhi);
            dim.n_phi = nPhi;            
        case 'X0'
            [options.priors.muX0,options.priors.SigmaX0, stateLabel] = feval(@priorUglyfier, varargin{:});
            options.inF.stateLabel =  stateLabel;
            options.inG.stateLabel =  stateLabel;
            nX0 = length(options.priors.muX0);
            dim.n = nX0;
    end