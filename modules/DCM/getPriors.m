function [priors] = getPriors(nreg,n_t,options,reduced_f,stochastic)
% builds priors for DCM inversion
% function [priors] = getPriors(nreg,n_t,options,reduced_f,stochastic)
% IN:
%   - nreg: # regions in the network
%   - n_t: length of the time series
%   - options: options structure, containing at least the .inF and .inG
%   fields, as built using prepare_fullDCM.m
%   - reduced_f: a flag for fixing a subpart of the hemodynamic parameters
%   (oxygen extraction fraction at rest, vasodilatory signal feedback rate
%   and vessel stifness) to their prior value
%   - stochastic: a flag for stochastic DCM
% OUT:
%   - priors: the 'priors' structure that can be used to invert the DCM
%   using VBA_NLStateSpaceModel.m

extended = isfield(options.inF,'extended') && options.inF.extended ;


%% get dimensions
try
    dim.n_theta = options.inF.ind5(end);
catch
    dim.n_theta = 0;
end

if extended
    dim.n_phi = options.inG.indr(end);
else
    dim.n_phi = options.inG.ind2(end);
end
dim.n = 5*nreg;


if extended
    n_r = numel(options.inG.r) ;
else
    n_r=0;
end

%% SET PRIORS

%% == initial conditions

nx = 5*nreg;
if extended
    nx=nx+n_r;
end
%- state %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO Different for predictors ?
priors.muX0 = zeros(nx,1); 
priors.SigmaX0 = 1e-4*eye(nx);


%% == evolution parameters
priors.muTheta = zeros(dim.n_theta,1); 

%- HRF
priors.SigmaTheta = 1e-2*eye(dim.n_theta);
if reduced_f
    % fix some HRF params to their default values
    priors.SigmaTheta(options.inF.ind1,options.inF.ind1) = 0;
    priors.SigmaTheta(options.inF.ind3,options.inF.ind3) = 0;
    priors.SigmaTheta(options.inF.ind5,options.inF.ind5) = 0;
end
%- DCM
idx = 1:(options.inF.indself-1);
priors.SigmaTheta(idx,idx) =   1e0*eye(length(idx));
priors.SigmaTheta(options.inF.indself,options.inF.indself) = 1e-1; %0.1
%- extension
if extended
    idx = options.inF.indself+1:options.inF.indhself-1;
    priors.SigmaTheta(idx,idx) =   1e0*eye(numel(idx)); %5
    priors.SigmaTheta(options.inF.indhself,options.inF.indhself) = 1e-1*eye(n_r) ; %.5
    % const
%     priors.muTheta(options.inF.indconst) = 0; 
%     priors.SigmaTheta(options.inF.indconst,options.inF.indconst) = 10*ones(dim.n_r); 
    
end


%% == observation parameters
%- HRF
priors.muPhi = zeros(dim.n_phi,1);
priors.SigmaPhi = 1e-2*eye(dim.n_phi);
%- extension
if extended
   priors.SigmaPhi(options.inG.indr,options.inG.indr)=20*eye(numel(options.inG.indr)) ; 
end


%% == state and measurement noise covariances
for t = 1:n_t
    dq = 1e2*ones(dim.n,1); 
    dq(options.inF.n5) = 1;
    if extended
        dq(options.inF.r) = 100;
    end
    priors.iQx{t} = diag(dq);         
end

% muxer
try
    gsi = find ([options.sources.type] == 0);
    n_Gsources = numel(gsi);
    for n=1:n_Gsources
        for t = 1 : n_t
            dim_y = numel (options.sources(gsi(n)).out);
            priors.iQy{t,n} = eye(dim_y);         
        end
    end
catch   
    n_Gsources = 1;
    for t = 1 : n_t
        priors.iQy{t,1} = eye(nreg);         
    end
end

%= precision hyperparameters
priors.a_sigma = 1e0*ones(n_Gsources,1);
priors.b_sigma = 1e0*ones(n_Gsources,1);
if ~stochastic
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
else
    TR = options.decim.*options.inF.deltat;
    priors.b_alpha = 1e0;%TR/1e2;
    priors.a_alpha = 1e0;
end


