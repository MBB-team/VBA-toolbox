
 %---------------
 
 load sources_betas

 
 L = [];
for IS = 1:32
 

Y = [];
U = [];

    Nsess = 6;
    Nsubject = 1;

    for i_subject = IS
for i_session = 1: Nsess


    y = Betas{i_subject}{i_session}.M1;
    u = [Betas{i_subject}{i_session}.leftTPJ;...
        Betas{i_subject}{i_session}.mpfc;...
        Betas{i_subject}{i_session}.caudate];

    U = [U;u];
    Y = [Y;y];
    
end
    end
 

 %[ f_fname_e,g_fname_e,dim_e,options_e ] = getModelforSessions( f_fname,g_fname,dim,options,in_sessions);
 
 %---------------


Ntrials = length(Y);


dim = struct('n',0,...  %
             'p',1,... % output (la combinaison)
             'n_theta',0,... % evolution parameters
             'n_phi', 3,... % observation parameters
             'n_t',Ntrials,...
             'u',3);



options.dim = dim;

g_fname = @g_lm;
f_fname = [];


in_sessions = struct();

in_sessions.ind.phi = ones(Nsess,1)*[1:3];
in_sessions.dim_e.n_phi = 3;
in_sessions.n_sess = 6;
in_sessions.f_fname = [];
in_sessions.g_fname = @g_lm;
in_sessions.inG = [];
in_sessions.inF = [];

in_sessions.binomial = 0;

[ f_fname_su,g_fname_su,dim_su,options_su ] = getModelforSessions( f_fname,g_fname,dim,options,in_sessions)


%%
% in_sessions = struct();
% in_sessions.dim_e.n_phi = 3*Nsubject;
% in_sessions.ind.phi = zeros(Nsess*Nsubject,3);
% for s = 1 : Nsubject
% in_sessions.ind.phi( ((s-1)*Nsess+1):s*Nsess ,:) = ones(Nsess,1)*[3*(s-1)+1:3*s];
% end 
in_sessions = struct();
in_sessions.dim_e.n_phi = 3;
in_sessions.ind.phi = zeros(Nsess*Nsubject,3);
for s = 1 : Nsubject
in_sessions.ind.phi( ((s-1)*Nsess+1):s*Nsess ,:) = ones(Nsess,1)*[1:3];
end 

in_sessions.n_sess = Nsubject;
in_sessions.f_fname = [];
in_sessions.g_fname = g_fname_su;
in_sessions.inG = options_su.inG;
in_sessions.inF = options_su.inG;
in_sessions.binomial = 0;

[ f_fname_pop,g_fname_pop,dim_pop,options_pop ] = getModelforSessions( f_fname_su,g_fname_su,dim_su,options_su,in_sessions)


%%
% -- declaring default priors
dim = dim_pop;
% Priors on parameters (mean and Covariance matrix)
priors.muPhi = zeros(dim.n_phi,1); 
priors.muTheta = zeros(dim.n_theta,1);
priors.SigmaPhi = 1e0*eye(dim.n_phi);
priors.SigmaTheta = 1e0*eye(dim.n_theta);
% No state noise for deterministic update rules
priors.a_alpha = Inf;
priors.b_alpha = 0;


options=options_pop;
options.DisplayWin = 1;
options.GnFigs = 0;
options.binomial = 0; % Dealing with binary data
options.priors = priors;
options.isYout = zeros(Nsess*Nsubject,Ntrials); % Excluding data points



% -- declaring modulatory inputs

mod = struct();
mod.indu = [1:3];
mod.names ={'mod1','mod2','mod3'};

% -- declaring modulation parameters


mod.phi.indp = [1:3;... % index of param in param vector
                1:3]; % index of modulating input

            
            
            
% -- Launching the inversions

%[posterior,out] = VBA_NLStateSpaceModel(Y,U,f_fname,g_fname,dim,options);
[inversions,inv_order,p_m,LogEv] = VBA_Inversion_modulation(Y,U,f_fname_pop,g_fname_pop,dim_pop,options,mod);

%[inversions,inv_order,p_m,LogEv] = VBA_Inversion_modulation(Y,U,f_fname,g_fname,dim,options,mod);

L = [L,LogEv];
end

%save raphael L inversions inv_order p_m

%%
group_level_analysis(L([6,8],:)','RFX') 

