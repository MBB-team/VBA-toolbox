function [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options,in)

% VB inversion of nonlinear stochastic DCMs
% function [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options,in)
% This function inverts any nonlinear state-space model of the form:
%   y_t = g( x_t,u_t,phi ) + e_t
%   x_t = f( x_t-1,u_t,theta ) + f_t
% using a variational bayesian annealing scheme.
% It is devoted to the trial estimation problem, ie:
%   - estimation of the dynamic hidden-states x_t
%   - estimation of the static observation/evolution parameters phi/theta
%   - estimation of the variance parameters associated with the measurement
%   error e_t and the stochastic innovations f_t.
%
% It is possible to specify priors that are such that the evolution
% function is deterministic. In such a case, the function inverts a
% standard DCM, and only evolution/observation parameters are estimated
% (see bellow).
%
% Calling this function without any input outputs the default inputs.
% NB: You can re-display the standard graphical output of the VB scheme by
% calling the VBA_ReDisplay.m function.
%
% IN :
%   - y: pxn_t mesurements matrix
%   - u: mxn_t known input matrix (which is required as an argument in
%   the obvservation/evolution functions, default is empty [])
%   - f_fname (resp. g_fname): name/handle of the function that outputs the
%   evolution (resp. observation) of the hidden states.
%       ! NB: The generic i/o form of these functions has to conform to the
%       following:
%                  [ fx,dfdx,dfdP,d2fxdP ] = fname(x_t,P,u_t,in)
%                      |________________|
%                         { varargout }
%       where:
%       . x_t, u_t and P are respectively the current hidden state value,
%       the current input value and the evolution/observation parameters
%       value
%       . in contains any additional information which has to be passed to
%       the functions. For example, it is used to indicate the chosen time
%       discretization (dt) for evolution functions that derive from a
%       continuous formulation.
%       . fx is the current evaluation of the function
%       . dfdx and dfdP are the matrices containing the gradients of the
%       function w.r.t. the hidden states and the invariant parameters
%       . d2fxdP is the nxn_phixp tensor containing the mixed derivatives.
%       ! NB2: the last dimension always refers to the dimension of the
%       vector-function, ie for the observation function:
%           dfdx(:,i) is the vector gradient of g(i) w.r.t. X
%           d2fdxdP(:,:,i) is the matrix of mixed derivatives of g(i)
%           w.r.t. X (columns) and Phi (rows):
%           d2fdxdP(:,:,i) =
%               [ d2G(i)/dX(1)dPhi(1) d2G(i)/dX(1)dPhi(2) ... d2G(i)/dX(1)dPhi(n_phi)
%                 d2G(i)/dX(2)dPhi(1) d2G(i)/dX(2)dPhi(2) ... d2G(i)/dX(2)dPhi(n_phi)
%                 ...
%                 d2G(i)/dX(n)dPhi(1) d2G(i)/dX(n)dPhi(2) ... d2G(i)/dX(n)dPhi(n_phi) ]
%       ! NB3: these mixed derivatives are unequal for non-continuous
%       functions (d2G_dXdPhi ~= d2G_dPhidX) !!!
%       ! NB4: the output of the evolution/observation functions does not
%       have to provide all above derivatives. For, e.g., the evolution
%       function, the possible output are: [fx], [fx,dfdx], [fx,dfdx,dfdP]
%       or [fx,dfdx,dfdP,d2fdxdP]. The missing derivatives are
%       automatically computed using numerical derivation.
%   - dim: a structure variable containing the dimensions of the 3 sets of
%   the model's unknown variables:
%       .n: the dimension of the hidden-states
%   	.n_theta: the dimension of the vector of evolution parameters theta
%   	.n_phi: the dimension of the vector of observation parameters phi
%       ! NB: the n_t and p dimensions are derived from the provided
%        observation matrix 'y' !
%   - options: user-defined structure containing specific informations
%   regarding the model, ie (see VBA_check.m):
%       .priors: a structure variable containing the priors sufficient
%       statistics over the evolution/observation/precision parameters and
%       hidden-states initial condition (see VBA_priors.m for standard
%       output)
%       ! NB: if the prior about the stochastic innovations precision is a
%        Dirac delta (priors.a_alpha = Inf, priors.b_alpha = 0), the model
%        becomes a deterministic (ODE) state-space model. This is used
%        during the initialization of the posterior.
%       ! NB2: time-dependent covariance structure for both the measurement
%       and the state noise can be passed to the inversion routine through
%       the .iQx and .iQy fields (these are cell aray of time-dependent
%       precision matrices).
%       .decim: the number of times the evolution function is applied
%       between each time sample {1}. This is used to increase the
%       micro-time resolution
%       .microU: a flag determining whether the input u is specified in
%       micro-resolution time (1) or in data sampling time ({0}). This is
%       useful for providing sub-sampling sparse inputs to the system.
%       .inF: a (possibly structure) variable containing the additional
%       (internal) fixed parameters which may have to be sent to the
%       evolution function {[]}
%       .inG: idem for the observation function {[]}
%       .checkGrads: a flag for eyeballing provided analytical gradients of
%       the evolution/observation function against automated numerical
%       derivatives
%       .updateX0: a flag indicating whether or not the initial conditions
%       should be updated ({1}:yes, 0:no).
%       .updateHP: a similar flag for the (precision) hyperparameters {1}
%       .backwardLag: a positive integer that defines the size of the
%       short-sighted backward pass for the hidden states VB upodate {1}.
%       NB: if 0, there is no backward pass.
%       .MaxIter: maximum number of VB iterations {32}
%       .MinIter: minimum number of VB iterations {1}
%       .TolFun: minimum absolute increase of the free energy {2e-2}
%       .DisplayWin: flag to display (1) or not (0) the main display
%       window {1}
%       .ignoreMF: binary variable allowing to account for (1) or discard
%       (0) mean-field additional terms in the update equations {1}
%       .gradF: binary variable that flags whether the regularized
%       Gauss_newton update scheme optimizes the free energy (1) or the
%       variational energy ({0}) of each mean-field partition.
%       .GnMaxIter: maximum number of inner Gauss-Newton iterations {32}
%       .GnTolFun: minimum relative increase of the variational energy for
%       the inner regularized Gauss-Newton loops {1e-5}
%       .GnFigs: flag to display (=1) or not (=0) the Gauss-Newton inner
%       loops display figures {0}.
%       .init0: flag for evolution/observation parameters posterior
%       initialization (1: use priors, {0}: use ODE fit).
%       .Laplace: flag for equilibrium free energy evaluation ({1})
%       .delays: dim.nX1 vector containing the discrete delays for the
%       evolution function.
%       NB: nonzero delays induces an embedding of the system, thereby
%       increasing the computational demand of the inversion. The state
%       space dimension if multiplied by max(options.delays)+1!
%       .binomial: {0}. If 1, it is assumed that the likelihood function is
%       a binomial pdf, whose probability is given by the observation
%       function, i.e. g(phi) = p(y=1|phi). The Laplace approximation on
%       observation parameters still holds and the i/o of the inversion
%       routine is conserved.
%       .figName: the name of the display window
%   - in: structure variable containing the output of a previously ran
%   VBA_NLStateSpaceModel.m routine. Providing this input to the function
%   allows one to re-start the VB updates from the point where it was
%   stopped. The required fields are:
%       .posterior: the posterior structure variable (see above)
%       .out: the additional output structure of this routine (see above)
%
% OUT:
%   - posterior: a structure variable whose fields contains the sufficient
%   statistics (typically first and second order moments) of the
%   variational approximations of the posterior pdfs over the
%   observation/evolution/precision parameters and hidden-states time
%   series. Its fields are:
%   	.muX: posterior mean of the hidden states X (nxn_t matrix)
%       .SigmaX: covariance matrices of the variational posterior pdf of
%       the dynamic hidden-states. Using the Kalman-Rauch marginalization
%       procedure, this is further divided into:
%           SigmaX.current{t} : the instantaneous covariance matrix at t
%           SigmaX.inter{t} : the lagged covariance matrix between instants
%           t and t+1
%       .muX0: posterior mean of the hidden-states initial condition, ie
%       before the first observation (nx1 vector)
%       .SigmaX0: covariance matrix of the Gaussian df over the
%       hidden-states initial condition (nxn matrix)
%   	.muTheta: posterior mean of the evolution parameters (vector)
%       .SigmaTheta: covariance matrix of the variational posterior pdf of
%       the static evolution parameters
%   	.muPhi: posterior mean of the observation parameters (vector)
%       .SigmaPhi: covariance matrix of the variational posterior pdf of
%       the static observation parameters
%       .a_alpha / .b_alpha: shape and scale parameters of the variational
%       posterior pdf of the stochastic innovations precision
%       .a_sigma / .b_sigma: shape and scale parameters of the variational
%       posterior pdf of the measurement noise precision
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


% JD, 26/02/2007

warning off
options.tStart = tic;

%------------------- Dummy call ------------------------%
if nargin == 0
    disp('Outputing the default inputs to the VBA_NLStateSpaceModel.m routine.')
    disp('These are stored in the ''out'' variable.')
    [y,u,f_fname,g_fname,dim,options,in] = VBA_getDefaults;
    posterior = [];
    out.y = y;
    out.u = u;
    out.f_fname = f_fname;
    out.g_fname = g_fname;
    out.dim = dim;
    out.options = options;
    out.in = in;
    disp(out)
    return
elseif isweird(y)
    disp('Error: there is a numerical trouble with provided data!')
    posterior = [];
    out = [];
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
        F = in.out.F;
        if ~options.OnLine
            if isequal(options.g_fname,@VBA_odeLim)
                VBA_disp('Initialization: non-stochastic system (ODE).',options)
            else
                VBA_disp('Skipping initialization.',options)
                VBA_disp('Main VB inversion...',options)
            end
        end
    catch
        disp('Error: the ''in'' structure was flawed...')
        posterior = [];
        out = [];
        return
    end
    
else

    % Check input arguments consistency
    [options,u,dim] = VBA_check(y,u,f_fname,g_fname,dim,options);
    VBA_disp(' ',options)
    % NB: priors are built using VBA_priors.m during the compilation
    % of the function VBA_check.m. These are checked and stored in
    % options.priors.

    % Initialize posterior pdf
    try
        [posterior,suffStat,options] = VBA_Initialize(y,u,f_fname,g_fname,dim,options);
    catch
        disp('')
        disp('Error: could not initialize VB scheme: check priors!')
        posterior = [];
        out = [];
        return
    end
        
    suffStat.u = u;
    % NB: when inverting a full state-space model, the initialization
    % actually inverts its deterministic variant, i.e. an ODE-like
    % state-space model.

    % Store free energy after the initialization step
    F = options.init.F;
    
    if ~options.OnLine
        VBA_disp('Main VB inversion...',options)
    end

end

suffStat.F = F;

% clear persistent variables
clear VBA_odeLim
clear VBA_smoothNLSS

% loop variables
stop = 0; % flag for exiting VB scheme
it = 0; % index of VB iterations

% Display initialized posterior
[options] = VBA_initDisplay(options);
try
    delete(intersect(findobj('tag','diagnostics_tabs'),get(options.display.hfp,'children')));
end
VBA_updateDisplay(F,posterior,suffStat,options,y,0,'precisions')
if dim.n_phi > 0
    VBA_updateDisplay(F,posterior,suffStat,options,y,0,'phi')
end
if dim.n > 0
    VBA_updateDisplay(F,posterior,suffStat,options,y,0,'X')
end
if dim.n_theta > 0
    VBA_updateDisplay(F,posterior,suffStat,options,y,0,'theta')
end

%------------------------------------------------------------%
%----------------- Main VB learning scheme ------------------%
%------------------------------------------------------------%

while ~stop
    
    
    it = it +1; % Iteration index

    %------ Gauss-Newton VB IEKS (+ initial conditions) ------%
    if dim.n > 0
        % Call for variational Kalman-Rauch smoothing algorithm.
        % This updates both posterior and suffStat variables
        [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,'X');
        % Call for VB initial condition update.
        if options.updateX0
            [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,'X0');
        end
    end

    %---------- Observation parameters ---------%
    if dim.n_phi > 0
        % Call for Laplace-VB update rule for the observation
        % parameters. This updates the posterior and suffStat variables
        [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,'Phi');
    end

    %----------- Evolution parameters ----------%
    if dim.n_theta > 0
        % Call for Laplace-VB update rule for the observation
        % parameters. This updates the posterior and suffStat variables
        [posterior,suffStat] = VBA_GN(y,posterior,suffStat,dim,u,options,'Theta');
    end

    %---------- Precision parameters ---------------%
    % Call update VB rule for precision parameters. This updates posterior
    % and suffStat variables.
    if options.updateHP
        [posterior,suffStat] = VBA_VBVarParam(y,posterior,suffStat,dim,options);
        [Fi] = VBA_FreeEnergy(posterior,suffStat,options);
        suffStat.F = [suffStat.F,Fi];
    else
        Fi = suffStat.F(end);
    end
    
    % Display progress
    VBA_updateDisplay(suffStat.F,posterior,suffStat,options,y,it,'precisions')
    VBA_updateDisplay([F Fi],posterior,suffStat,options,y,it,'F')
    
    %--------------- Termination condition ---------------%
    dF = Fi - F;
    F = Fi;
    if (  ( (abs(dF)<=options.TolFun)||it==options.MaxIter ) &&  it >=options.MinIter  ) || isweird(dF)
        stop  = 1;
        [Fi] = VBA_FreeEnergy(posterior,suffStat,options);
        suffStat.F = [suffStat.F,Fi];
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

[posterior,out] = VBA_wrapup(posterior,options,dim,suffStat,u,y,it);


