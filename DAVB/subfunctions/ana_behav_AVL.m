function [posterior,out] = ana_behav_AVL(y,u,flag,iQy)
% OTO: inverts audio-visual associative learning model given reaction time data
% function [posterior,out] = ana_behav_AVL(y,u,flag,iQy)
% IN:
%   - y: 1xt row vector of reaction times
%   - u: 2xt matrix of (visual) outcome and subjects categorical decisions.
%   NB: u(1,:) should be 1 when outcome='house' and -1 when outcome='face'
%   and u(2,:) should be 1 when choice='house' and 0 when choice='face'
%   - flag: perceptual model (1=static, 2=dynamic, 3=volatile)
%   - iQy: 1xt cell array of the data precision. This can be used to
%   discard tirls (e.g. perceptual categorization errors).
% OUT:
%   - posterior: structure containing the conditional sufficient statistics
%   of the audio-visual asociative learning model parameters
%   - out: model inversion infos


if ~exist('flag','var') || isempty(flag)
    flag = 2;
end
if ~exist('iQy','var') || isempty(iQy)
    n_t = size(y,2);
    priors.iQy = cell(n_t,1);
    for t= 1:n_t
        priors.iQy{t} = 1;
    end
else
    priors.iQy = iQy;
end


inF.n = 8;          % max # iterations (per trial) for the VB observer
inF.uu = 1;         % index of sensory signals in the vector u
inG.uc = 2;         % index of observer's choices in the vector u
f_fname = @f_AVL;
g_fname = @g_AVL;
switch flag
    case 1
        priors.muX0 = [0.5;0;1e3;0];
        priors.SigmaX0 = 0*eye(4);
        priors.SigmaX0(3,3) = 1e2;
        priors.muTheta = [2;0];
        priors.SigmaTheta = 0*1e-2*eye(2);
    case 2
        priors.muX0 = [0.5;0;1e3;0];
        priors.SigmaX0 = 0*eye(4);
        priors.muTheta = [2;0];
        priors.SigmaTheta = 1e2*eye(2);
        priors.SigmaTheta(1,1) = 0;
    case 3
        priors.muX0 = [0.5;0;1e0;-2;1e0;0];
        priors.SigmaX0 = 0*eye(6);
        priors.muTheta = [0;-2];
        priors.SigmaTheta = 1e1*eye(2);
        priors.SigmaTheta(1,1) = 0;
end
priors.muPhi = [0;0];%[1;-1];
priors.SigmaPhi = 1e2*eye(length(priors.muPhi));
priors.a_alpha = Inf;
priors.b_alpha = 0;
priors.a_sigma = 1e2;
priors.b_sigma = 1e-2;

% Build options and dim structures for model inversion
options.priors      = priors;
options.inF         = inF;
options.inG         = inG;
options.GnFigs      = 0;
options.DisplayWin  = 0;
options.inF.flag    = flag;
% options.gradF       = 1;
dim.n_phi           = length(priors.muPhi);
dim.n_theta         = 2;
switch options.inF.flag
    case 1
        dim.n       = 4;
    case 2
        dim.n       = 4;
    case 3
        dim.n       = 6;
end

% call OTO model inversion routine
[posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);




