function [theta,phi] = get_HRFparams(TR,microDT,lin,verbose)

% get Ballon model parameters from canonical SPM HRF

% function [theta,phi] = get_HRFparams(TR,microDT,lin,verbose)
% IN:
%   - TR: the data sampling (in sec)
%   - microDT: the micro-time resolution (in sec)
%   - lin: a flag for the linearized version of the Balloon model {0}
%   - verbose: a flag for displaying (1) or not ({0}) the inversion results
% OUT:
%   - theta: fitted evolution parameters
%   - phi: fitted observation parameters

% try, [theta,phi] = defaultHRFparams; return; end

if nargin <3
    lin = 0;
else
    lin = ~~lin;
end
if nargin <4
    verbose = 0;
else
    verbose = ~~verbose;
end

%-- First initialize with finer sampling frequency if required
if TR > 1
    [muTheta] = get_HRFparams(1,microDT,lin,verbose);
end

%-- Then invert HRF model
% get spm canonical hrf
[hrf] = VBA_spm_hrf(TR);
% get basic i/o
n_t         = length(hrf);          % number of time samples (over 40 sec)
decim       = max([1,round(TR./microDT)]);
u           = zeros(1,n_t*decim);         % input
u(1,1)      = 1e0;
f_fname     = @f_HRF2;               % Ballon model evolution function
g_fname     = @g_HRF3;               % Ballon model observation function
inF.linearized = lin;
inF.deltat  = microDT;
% Build priors for model inversion
priors.muX0         = [0;0;0;0];
priors.SigmaX0      = 0e-2*eye(4);
try
    % if initialized with TR = 1
    priors.muTheta = muTheta;
catch
    priors.muTheta      = 0*ones(6,1);
end
priors.SigmaTheta   = 1e-3*eye(6,6);
priors.SigmaTheta(5,5) = 0;
priors.muPhi        = 0*ones(2,1);
priors.SigmaPhi     = 1e-3*eye(2,2);
priors.a_alpha      = Inf;
priors.b_alpha      = 0;
priors.a_sigma      = 1e6;
priors.b_sigma      = 1e0;
% Build options and dim structures for model inversion
options.priors      = priors;
options.DisplayWin  = verbose;
options.verbose     = verbose;
options.inF         = inF;
options.decim       = decim;
options.microU      = 1;
options.inG         = [];
dim.n_theta         = 6;
dim.n_phi           = 2;
dim.n               = 4;
% call inversion routine
[posterior] = VBA_NLStateSpaceModel(hrf',u,f_fname,g_fname,dim,options);

% refine undershoot fit
% priors.muTheta = posterior.muTheta;
% priors.muPhi = posterior.muPhi;
% for t=1:n_t
%     if (t>=15 && t<30)
%         priors.iQy{t} = 1e1;
%     elseif (t>=4 && t<8)
%         priors.iQy{t} = 1e1;
%     else
%         priors.iQy{t} = 1e0;
%     end
% end
% options.priors = priors;
% [posterior] = VBA_NLStateSpaceModel(hrf',u,f_fname,g_fname,dim,options);

% extract Ballon model evolution/observation parameters
theta = posterior.muTheta;
phi = posterior.muPhi;

% if verbose
%     disp('!Extracting canonical HRF parameters: pause!')
%     pause
% end



