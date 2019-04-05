function [posterior,out] = VBA_UNL0(y,u,g_fname,dim,options,in)

% VB inversion of generative models with unnormalized likelihoods
% function [posterior,out] = VBA_UNL0(y,u,g_fname,dim,options,in)
% This function inherits its i/o format from the main VBA model inversion
% routine VBA_NLStateSpaceModel.m, but focuses on static models (no hidden
% states nor evolution function).
% IN :
%   - y: 1Xn_t mesurements vector
%   - u: mXn_t known input vector
%   - g_fname: name/handle of the function that returns the unnormalized
%   likelihood evaluated at the data y.
%   - dim: a structure (for i/o compatibility with VBA_NLStateSpaceModel)
%   containing the dimensions of the set of unknown model variables:
%   	.n_phi: the dimension of the vector of observation parameters phi
%   - options: user-defined structure containing specific informations
%   regarding the model, ie (see VBA_NLStateSpaceModel)
%   - in: structure variable containing the output of a previously ran
%   VBA_NLStateSpaceModel.m routine. Providing this input to the function
%   allows one to re-start the VB updates from the point where it was
%   stopped. The required fields are:
%       .posterior: the posterior structure variable (see above)
%       .out: the additional output structure of this routine (see above)
%
% OUT:
%   - posterior: a structure variable (for i/o compatibility with
%   VBA_NLStateSpaceModel) whose fields contain the sufficient statistics
%   of the variational approximations to the posterior pdfs over the
%   observation parameters:
%   	.muPhi: posterior mean of the observation parameters (vector)
%       .SigmaPhi: covariance matrix of the variational posterior pdf of
%       the static observation parameters
%   - out: a structure variable containing the fields...
%       .CV: convergence flag (0 if the algorithm has stopped because it
%       reached the options.MaxIter termination condition)
%       .F: the free energy associated with the inversion of the model
%       .options: the options structure which has been used for the
%       model inversion
%       .dim: the dimensions of the model
%       .it: the number of iterations which have been required for reaching
%       the convergence criteria
%       .suffStat: a structure containing internal variables that act as
%       sufficient statistics for the VB updates (e.g. predicted data ...)
%       .fit: structure, containing the following fields:
%           .LL: log-likelihood of the model
%           .R2: coefficient of determination (the fraction of variance
%           unexplained is 1-R2)
%           .AIC: Akaike Information Criterion
%           .BIC: Bayesian Informaion Criterion


% JD, 09/03/2016

warning off
options.tStart = tic;
options.UNL = 1;
posterior = [];
out = [];

%------------------- Dummy call ------------------------%
if VBA_isWeird(y)
    disp('Error: there is a numerical trouble with provided data!')
    return
end

%-------------------- Initialization -------------------%

if exist('in','var')
    
    try  % Skip checking and replace prior initialization
        dim = in.out.dim;
        posterior = in.posterior;
        suffStat = in.out.suffStat;
        options = in.out.options;
        u = in.out.u;
        if ~options.OnLine
            VBA_disp('Skipping initialization.',options)
            VBA_disp('Main VB inversion...',options)
        end
    catch
        disp('Error: the ''in'' structure was flawed...')
        return
    end
    
else
    
    % Check input arguments consistency (and fill in priors if necessary)
    % -- !! This forces static model inversion !! --
    f_fname = @f_Id; % dummy evolution function (for VBA_check)
    dim.n = 0;
    dim.n_theta = 0;
    dim.n_t = size(y,2);
    dim.p = 1;
    % -- --
    [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options);
    VBA_disp(' ',options)
    
    try % Initialize posterior pdf
        posterior = options.priors;
        [Iphi,SigmaPhi,deltaMuPhi,suffStat] = VBA_Iphi_UNL(posterior.muPhi,y,posterior,[],dim,u,options);
        options.init.F = VBA_FreeEnergy_UNL(posterior,suffStat,options);
        options.init.posterior = posterior;
        options.init.suffStat = suffStat;
    catch
        VBA_disp(' ',options)
        VBA_disp('Error: Program stopped during initialisation',options)
        return
    end
    
    % Store free energy after the initialization step
    suffStat.F = options.init.F;
    
    VBA_disp('Main VB inversion...',options)
    
end



% Display initialized posterior
[options] = VBA_initDisplay(options);
% delete(options.display.ha(6))
try
    delete(intersect(findobj('tag','diagnostics_tabs'),get(options.display.hfp,'children')));
end
if dim.n_phi > 0
    VBA_updateDisplay(posterior,suffStat,options,y,0,'phi')
end
if options.updateHP
    VBA_updateDisplay(posterior,suffStat,options,y,0,'precisions')
end

%------------------------------------------------------------%
%----------------- Main VB learning scheme ------------------%
%------------------------------------------------------------%

it = 0; % index of VB iterations
stop = it>=options.MaxIter; % flag for exiting VB scheme

while ~stop
    
    it = it +1; % iteration index
    F0 = suffStat.F(end);
    
    %---------- Observation parameters ---------%
    if dim.n_phi > 0
        % Call for Laplace-VB update rule for the observation
        % parameters. This updates the posterior and suffStat variables
        [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,'Phi');
    end
    
    %---------- Temperature parameter ---------------%
    % Call update VB rule for temperature. This updates posterior
    % and suffStat variables.
    if options.updateHP
        [posterior,suffStat] = VBA_UNLtemp(y,posterior,suffStat,dim,u,options);
        suffStat.F = [suffStat.F,VBA_FreeEnergy(posterior,suffStat,options)];
        VBA_updateDisplay(posterior,suffStat,options,y,it,'precisions')
    end
    
    
    % Display progress
    if ~options.OnLine && options.verbose
        dF = suffStat.F(end)-F0; % FE increase during the last VB iteration
        fprintf(['VB iteration #',num2str(it),'         F=','%4.3e','         ... dF=','%4.3e'],suffStat.F(end),dF(end))
        fprintf('\n')
    end
    
    %--------------- Termination condition ---------------%
    dF = diff(suffStat.F);
    dF = dF(end);
    if (  ( (abs(dF)<=options.TolFun)||it==options.MaxIter ) &&  it >=options.MinIter  ) || VBA_isWeird(dF)
        stop  = 1;
        if abs(dF) <= options.TolFun
            out.CV = 1;
        else % did not converge
            out.CV = 0;
        end
    end
    
end



%------------------------------------------------------------%
%---------------- wrap up inversion results -----------------%
%------------------------------------------------------------%

suffStat.dx = [];
suffStat = VBA_Hpost(posterior,suffStat,options);
[posterior,out] = VBA_wrapup(posterior,options,dim,suffStat,u,y,it);


