function [y,u,f_fname,g_fname,dim,options,in]    = VBA_getDefaults()
% outputs default input to VBA_NLStateSpaceModel.m
% function [y,u,f_fname,g_fname,dim,options,in]    = VBA_getDefaults()
%
% This function provides the defaults inputs to the main
% VBA_NLStateSpaceModel.m routine. See VBA_check.m routine.

%------- Main inputs ------------%
y                   = 'pxn_t measured data matrix';
u                   = 'n_uxn_t driving input matrix (can be left empty)';
f_fname             = 'name/handle of the evolution function';
g_fname             = 'name/handle of the observation function';
dim.n               = 'dimension of the state space';
dim.n_theta         = 'dimension of the evolution parameter vector';
dim.n_phi           = 'dimension of the observation parameter vector';
in                  = 'only for inversion refinement (can be unspecified)';


%------- Options structure -------%
options.priors = 'structure containing the first two moments of all prior pdfs [see VBA_defaultPriors()]';
options.decim = 'for micro-time resolution';
options.microU = 'flag for micro-time input u';
options.inF = '(internal) parameters of the evolution function';
options.inG = '(internal) parameters of the observation function';
options.checkGrads = 'Check analytical gradients against numerical gradients';
options.updateX0 = 'Flag for initial conditions update';
options.updateHP = 'Flag for precision hyperparameters update';
options.MaxIter = 'Maximum number of iterations';
options.MinIter = 'Minimum number of iterations';
options.TolFun = 'Minimum change in the free energy ';
options.DisplayWin = 'VB display window';
options.GnMaxIter = 'Gauss-Newton inner-loops maximum number of iterations';
options.GnTolFun = 'Gauss-Newton inner-loops threshold on relative variational energy';
options.GnFigs = 'Gauss-Newton inner-loops figures';
options.verbose = 'summary messages in MATLAB main window';
options.finalEval = 'Externally-specified function evaluation: end of VB algo';
options.figName = 'Name of the display figure';
options.delays = 'dim.nX1 vector for delay embedding';
options.sources = 'Type and dimension of the data distribution(s)';



