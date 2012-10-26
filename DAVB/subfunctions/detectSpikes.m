function [posterior,out] = detectSpikes(Y,flag,nmax)

% analyzes calcium imaging datafile
% function [posterior,out] = detectSpikes(Y,flag,nmax)
% IN:
%   - Y: data structure
%   - flag: model type ('HH','fitzhugh','lognormal' or 'bumps')
%   - nmax: [only for 'bumps'] max nmber of input bump functions
% OUT:
%   - posterior,out: the output of the model inversion

if nargin < 2
    flag = 'bumps';
end
if nargin < 3
    nmax = 12;
end




% common
n_t = size(Y.calcium,2);
g_fname = @g_Id;
dt = mean(diff(Y.calcium_time));
2e-2/dt
options.decim = round(2e-2/dt);
y = Y.calcium(:,1:options.decim:n_t);
options.inF.delta_t = dt;
options.inF.sc              = 5e1; % time scale change for HH model
options.inF.a               = 0.5; % prior calcium temporal decay rate
options.inG.ind             = 1; % index of observable state

switch flag

    case 'HH'
        
        % this model uses Hodgkin-Huxley AP propagation model, convolved
        % with a 1st-order convolution kernel (calcium response).

        f_fname = @f_HH_calcium;
        g_fname = @g_Id;
        
        % Build priors for model inversion
        priors.muX0 = [0;0;-3;-3/4;1/3];
        priors.SigmaX0 = 0e0*eye(5);
        priors.muTheta = [1;1;0*ones(3,1)];
        priors.SigmaTheta = 0e1*eye(5);
        priors.SigmaTheta(1,1) = 1e0;
        priors.a_alpha      = 1e0;
        priors.b_alpha      = 1e0;
        priors.a_sigma      = 1e0;
        priors.b_sigma      = 1e0;
        for t = 1:n_t
            dq              = 1e2*ones(5,1);
            dq(2)           = 1e0;
            priors.iQx{t}   = diag(dq);
        end
        options.priors = priors;
        options.updateHP    = 1;
        options.backwardLag = 32;
        options.GnFigs = 0;
        options.decim = 1;
        options.MaxIterInit = 0;
        
        dim.n_theta         = 5;
        dim.n_phi           = 0;
        dim.n               = 5;
    
        u = [];
        
    
    case 'fitzhugh'
        
        % this model uses Fitz-Hugh-Nagumo AP propagation model, convolved
        % with a 1st-order convolution kernel (calcium response).
        
        f_fname = @f_FitzHughNagumo_calcium;
        g_fname = @g_Id;
    
        priors.muX0 = [0;-.2;0];
        priors.SigmaX0 = 1e-1*eye(3);
        priors.muTheta = 0.*ones(5,1);
        priors.SigmaTheta = 1e-3*eye(5);
        priors.SigmaTheta(2,2) = 0;
        priors.a_alpha      = 1e0;
        priors.b_alpha      = 1e0;
        priors.a_sigma      = 1e0;
        priors.b_sigma      = 1e0;
        for t = 1:n_t
            dq              = 1e2*ones(3,1);
            dq(2)           = 1;
            priors.iQx{t}   = diag(dq);
        end
        
        options.updateHP    = 0;
        options.backwardLag = 8;
        options.priors = priors;
        dim.n_theta         = 5;
        dim.n_phi           = 0;
        dim.n               = 3;
        
        u = [];

    
    case 'lognormal'
        
        % this model uses an exponentiated AR(1) model for AP-induced
        % membrane depolarization, convolved with a 1st-order convolution
        % kernel (calcium response).
        
        f_fname             = @f_calcium;
        
        % Build priors for model inversion
        priors.muX0         = [0;-.2];
        priors.SigmaX0      = 1e0*eye(2);
        priors.muTheta      = 0*ones(2,1);
        priors.SigmaTheta   = 1e0*eye(2);
        priors.a_alpha      = 1e0;
        priors.b_alpha      = 1e0;
        priors.a_sigma      = 1e0;
        priors.b_sigma      = 1e0;
        for t = 1:n_t
            dq              = 1e2*ones(2,1);
            dq(2)           = 1;
            priors.iQx{t}   = diag(dq);
        end
        
        % Build options and dim structures for model inversion
%         options.inF.a       = 1;
        options.decim       = 1;
        options.GnFigs      = 0;
        options.DisplayWin  = 1;
        options.MinIter     = 1;
        options.updateHP    = 1;
        options.backwardLag = 8;
        options.MaxIterInit = 0;
        options.GnMaxIter   = 16;
        dim.n_theta         = 2;
        dim.n_phi           = 0;
        dim.n               = 2;
        
        % input
        u                   = [];
        
    case 'bumps'
        
        f_fname             = @f_lin1D;
        
%         options.inF.a       = 1;
        options.inF.u_fname = @u_GaussianBumps;
        options.inF.indscale = 2:nmax+1;
        options.inF.indcentres = nmax+2:2*nmax+1;
        
        % Build priors for model inversion
        priors.muX0         = 0;
        priors.SigmaX0      = 0e-0*eye(1);
        priors.muTheta      = 0*ones(2+2*nmax,1);
        priors.muTheta(2)   = log(options.inF.delta_t)./4; % bumps log spread
        priors.muTheta(options.inF.indcentres) = ...
            log(options.inF.delta_t+...
            [0:nmax-1]*options.inF.delta_t*(n_t-1)./(nmax-1));
        priors.SigmaTheta   = 1e2*eye(2+2*nmax);
        priors.a_alpha      = Inf;
        priors.b_alpha      = 0;
        priors.a_sigma      = 1e0;
        priors.b_sigma      = 1e0;
        
        % Build options and dim structures for model inversion
        options.decim       = 1;
        options.GnFigs      = 0;
        dim.n_theta         = 2+2*nmax;
        dim.n_phi           = 0;
        dim.n               = 1;
        
        % input
        u                   = [zeros(1,n_t);...
            0:options.inF.delta_t:options.inF.delta_t*(n_t-1)];

    otherwise
        
        disp('unknown model')
        posterior = [];
        out = [];
        return
        
end


options.GnFigs = 1;
options.priors              = priors;
[posterior,out]             = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
VBA_ReDisplay(posterior,out);

