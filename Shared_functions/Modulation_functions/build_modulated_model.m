function [ f_fname_m,g_fname_m,dim_m,options_m, mod_m ] = build_modulated_model( f_fname,g_fname,dim,options,mod )
% Builds an extended model that is the modulated model.
% In the modulated model, parameters depend on inputs through
% linear combination.
% Parameters of the modulated model can be separeted in two sets
% - initial parameters p_s of the initial model
% - modulation parameters p_m
% A modulation parameter is associated to a modulation variable u.
% A modulation parameter is always associated to an initial parameter
% An initial parameter can be associated to many modulation parameters
% I call extended parameter p_e a group consisting of an initial parameter
% and its associated modulation parameters
% p_e = p_s + ( p_m(1)*u(1)+...p_m(K)*u(K) )
% In the abscence of influence of u on p_s, p_m =0;
%
% INPUT
% - f_fname : handle of the evolution function of the initial model
% - g_fname : handle of the observation function of the initial model
% - dim : dimension of the initial model
% - options : options of the initial model
% - mod : structure containing information about desired modulation
%   - .phi/theta : cell of initial parameters that are modulated
%     - .indu : indices of the modulatory inputs
% OUTPUTS
% - f_fname_m : handle of the evolution function of the modulated model
% - g_fname_m : handle of the observation function of the modulated model
% - dim_m : dimension of the modulated model
% - options : options of the modulated model
% - mod_m : structure containing information the modulation
%   - .phi/theta : cell of initial parameters that are modulated
%     - .indu : indices of the modulatory inputs
%     - .indpi : index of new initial parameter in the modulated model
%     - .indpm : index of mudulating parameters associted to inputs indu in
%     the modulated model

dim_m = dim;
f_fname_m = @f_mod;
g_fname_m = @g_mod;

mod_m = mod;

i_pm = 1; % index in modulated model
for i_p = 1 : dim.n_phi
    mod_m.phi{i_p}.indpi = i_pm;
    mod_m.phi{i_p}.indp = i_p;
    mod_m.phi{i_p}.indpm = [];

    i_pm = i_pm + 1;
    try
        dim_m.n_phi = dim_m.n_phi + length( mod.phi{i_p}.indu ); % update parameter dimensions
        indpm = [];
        for i_u = mod.phi{i_p}.indu
            indpm =[indpm,i_pm];
            i_pm = i_pm + 1;

        end
        mod_m.phi{i_p}.indpm = indpm;
    catch
    end
end


%------------------------------------------
% Evolution parameters

i_pm = 1; % index in modulated model
for i_p = 1 : length(mod_m.theta)
    mod_m.theta{i_p}.indpi = i_pm;
    mod_m.theta{i_p}.indp = i_p;
    mod_m.theta{i_p}.indpm = [];

    i_pm = i_pm + 1;
    try
        dim_m.n_theta = dim_m.n_theta + length( mod.theta{i_p}.indu ); % update parameter dimensions
        indpm = [];
        for i_u = mod.theta{i_p}.indu
            indpm =[indpm,i_pm];
            i_pm = i_pm + 1;

        end
        mod_m.theta{i_p}.indpm = indpm;
    catch
    end
end

%------------------------------------------


options_m = options;
options_m.dim = dim_m;

options_m.inF.dim = dim; % initial model
options_m.inG.dim = dim; % initial model

options_m.inF.f_fname = f_fname; % initial model
options_m.inG.g_fname = g_fname; % initial model

options_m.inF.inF = options.inF; % initial model
options_m.inG.inG = options.inG; % initial model

options_m.inF.theta = mod_m.theta; %
options_m.inG.phi = mod_m.phi; %
%
end





