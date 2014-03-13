function [] = demo_gradf()
% [broken]

warning on
warning('This requires a .mat file containing simulations parameters!')
return

clear variables
close all

% Choose basic settings for simulations
f_fname = @f_gradf;
g_fname = @g_Id;

lo = load('Fs.mat');

n = size(lo.F,1);

Theta = zeros(3,n);

for i=1:n

    i
    
    y = lo.F(i,:);
    
    priors.muX0 = 2*y(1);
    priors.SigmaX0 = 1e8;
    priors.muTheta = [0;y(end)];
    priors.SigmaTheta = speye(2);
    priors.SigmaTheta(2,2) = 0;
    priors.a_alpha = Inf;
    priors.b_alpha = 0;
    priors.a_sigma = 1e2;
    priors.b_sigma = 1e-2;
    
    
    dim.n = 1;
    dim.n_theta = 2;
    dim.n_phi = 0;
    options.priors = priors;
    options.GnFigs = 0;
    options.inF = [];
    options.inG = [];
    options.DisplayWin = 0;
    
    [posterior,out] = VBA_NLStateSpaceModel(y,[],f_fname,g_fname,dim,options);
    
    
    Theta(:,i) = [posterior.muTheta;posterior.muX0];

end
