function [ f_fname_m,g_fname_m,dim_m,options_m, mod_m ] = build_modulated_model( f_fname,g_fname,dim,options,mod )
% build the modulated model.
% modulated model has more parameters corresponding to the parameters
% modulating the parameters of the initial model

dim_m = dim;
f_fname_m = @f_mod;
g_fname_m = @g_mod;

mod_m = mod;

%------------------------------------------
% Observation parameters
i_pm = 1;
for i_p = 1 : length(mod_m.phi)
    mod_m.phi{i_p}.indpe = i_pm;
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

i_pm = 1;
for i_p = 1 : length(mod_m.theta)
    mod_m.theta{i_p}.indpe = i_pm;
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
end

