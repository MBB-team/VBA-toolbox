function [ f_fname_e,g_fname_e,dim_e,options_e,fb_e ] = makeExtendedModel(dim,options,in_sessions)
% This function creates the description files of an Extended Model (EM)
% of an Initial Model (IM).
% The output of the EM is the concatenation of independant outputs of several IM models
% which may share parameters (evolution or observation parameters)

% The purpose of this function in the context of behavioral experiments is
% to have a generative modei for multiple independant sessions.
% IM : generative models for a single session
% EM : generative model for multiple sessions


% INPUT
% - dim : dimensions of the IM
%    - .u : dimension of the data vector used
% - options : options of the IM
% - in_sessions : information about sessions
%       - .f_fname : evolution function of the IM
%       - .g_fname : observation function of the IM model
%       - .dim : dimensions of the variables of the model for a single
%       session (base model)
%       - .ind.theta : indices of the variable theta used for each session
%       - .ind.phi : indices of the variable phi used for each session



n_sess = in_sessions.n_sess;
%---------------------------------------------------
%-- Dimensions of the extended model

%-- Check existence of indices

try in_sessions.ind.theta;
catch err
      msg = sprintf('You need to provide input in_sessions.ind.theta so we know which evolution parameters are used for which session');
      error(msg);
end  
try in_sessions.ind.phi;
catch err
      msg = sprintf('You need to provide input in_sessions.ind.phi so we know which observation parameters are used for which session');
      error(msg);
end  


%------------------------------------------------------
%-- Filling session structure

% --- General information for all sessions
in = struct();
in.nsess = in_sessions.n_sess;

% --- Specific information for each session
i_gx = 0; % cumulated index of output
i_x = 0; % cumulated index of hidden states
i_u = 0; % cumulated index of inputs

for i = 1 : n_sess
    
    try dim_s = dim{i};
    catch; dim_s = dim;end
     in.sess(i).dim = dim_s;
    
    % Information about the evolution/obsevation function for each session
    try in.sess(i).f_fname = in_sessions.f_fname{i};
    catch; in.sess(i).f_fname = in_sessions.f_fname;end
    try  in.sess(i).g_fname = in_sessions.g_fname{i};
    catch ;in.sess(i).g_fname = in_sessions.g_fname;end
    
    
%     % Information about the evolution/obsevation function for each session
%     in.sess(i).f_fname = in_sessions.f_fname; % the function to be used for each session
%     in.sess(i).g_fname = in_sessions.g_fname;
%     
    % Information about indices of parameters, hidden states and output used by
    % each session
    in.sess(i).ind.x = i_x+1:i_x+dim_s.n;   i_x = i_x + dim_s.n;
    in.sess(i).ind.gx = i_gx+1:i_gx+dim_s.p;   i_gx = i_gx + dim_s.p;
    in.sess(i).ind.u = i_u+1:i_u+dim_s.u;   i_u = i_u + dim_s.u;
    
    % Information about the evolution/obsevation paramaters for each session
    if isempty(in_sessions.ind.theta)
        in.sess(i).ind.theta = [];
    else in.sess(i).ind.theta = in_sessions.ind.theta(i,:); end
    if isempty(in_sessions.ind.phi)
        in.sess(i).ind.phi = [];
    else in.sess(i).ind.phi = in_sessions.ind.phi(i,:); end
    
    % Information about the evolution/obsevation extra input for each session
    try in.sess(i).inG = in_sessions.inG{i};
    catch; in.sess(i).inG = in_sessions.inG;end
    try  in.sess(i).inF = in_sessions.inF{i};
    catch ;in.sess(i).inF = in_sessions.inF;end
    
end

%------------------------------------------------------
%--- compute dimensions of the extended model

dim_e = struct();
dim_e.n_sess = n_sess;
dim_e.n_t = dim.n_t; % number of trials (unchanged)

 dim_e.p = 0; % output
 dim_e.u = 0; % number of trials (unchanged)
 dim_e.n = 0; % hidden states
 dim_e.n_theta = max(max(in_sessions.ind.theta));
 dim_e.n_phi = max(max(in_sessions.ind.phi));

for i = 1 : n_sess
dim_e.p = max([dim_e.p; in.sess(i).ind.gx(:)]);
dim_e.u = max([dim_e.u; in.sess(i).ind.u(:)]);
dim_e.n = max([dim_e.n; in.sess(i).ind.x(:)]);
end

%---- Options
options_e = options; % copy all options then modify

options_e.inF = in;
options_e.inG = in;

options_e.inG.dim = dim_e; % requested to know size of output when generating it.
options_e.dim = dim_e;


% options_e.GnFigs = 0;
% try options_e.binomial = in_sessions.binomial;
% catch; options_e.binomial = 0; end % default is continuous data
% try  options_e.DisplayWin = in_sessions.DisplayWin;
% catch ; options_e.DisplayWin = 0; end

options_e.isYout = zeros(dim_e.p,dim_e.n_t); 


%---- Function handles

f_fname_e = @f_nsess;
g_fname_e = @g_nsess;

%
%---- Handle for simulation
fb_e = [];
try
    
    %fb = options.fb;
    fb_e = struct('h_fname',@h_nsess,...
        'nsess',n_sess);
    fb_e.inH.indy = [];
    fb_e.inH.indfb = [];
    i_fb = 0;
    for i = 1 : n_sess
        
        % possible different models for each session
        try dim_s = dim{i};
        catch; dim_s = dim;end
        try fb_s = in_sessions.fb{i};
        catch; fb_s = in_sessions.fb;end
        
        fb_e.inH.sess(i).h_fname = fb_s.h_fname; % handle corresponding to session
        fb_e.inH.sess(i).indy = in.sess(i).ind.gx(fb_s.indy);
        fb_e.inH.sess(i).indfb = i_fb + [1:length(fb_s.indfb)]'; % index in the extended feedback vector (not in u!)
        
        try   fb_e.inH.sess(i).inH =  in_sessions.inH{i}; % case feedbacks not concatenated
        catch;  fb_e.inH.sess(i).inH =  in_sessions.inH(fb_e.inH.sess(i).indfb,:); end % case feedbacks concatenated
        
        
        fb_e.inH.indy = [fb_e.inH.indy; in.sess(i).ind.u(fb_s.indy)];
        fb_e.inH.indfb = [fb_e.inH.indfb; in.sess(i).ind.u(fb_s.indfb)];
        i_fb = i_fb+length(fb_s.indfb);

    end
    
    fb_e.inH.nsess = n_sess; % number of sessions
    fb_e.indy =fb_e.inH.indy;
    fb_e.indfb =fb_e.inH.indfb;
    
    
catch;
end
    options_e.fb = fb_e;

end






