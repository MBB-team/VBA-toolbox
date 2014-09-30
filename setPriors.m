function options = setPriors(options,priorType,varargin)

    switch priorType
        case 'theta'
            [options.priors.muTheta,options.priors.SigmaTheta, options.inF.paramLabel] = feval(@priorUglyfier, varargin{:});        
        case 'phi'
            [options.priors.muPhi,options.priors.SigmaPhi, options.inG.paramLabel] =  feval(@priorUglyfier, varargin{:});
        case 'X0'
            [options.priors.muX0,options.priors.SigmaX0, stateLabel] = feval(@priorUglyfier, varargin{:});
            options.inF.stateLabel =  stateLabel;
            options.inG.stateLabel =  stateLabel;
    end