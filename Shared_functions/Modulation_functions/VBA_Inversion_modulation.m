function [inversions,mod,LogEv] = VBA_Inversion_modulation(y,u,f_fname,g_fname,dim,options,mod,in)
% This function performs model comparision over a set of models.
% Models considered only differ in their priors, they share the same structure.
% Parameters that can be separeted in two sets
% - standard parameters p_s
% - modulation parameters p_m
% A modulation parameter is associated to a modulation variable u.
% A modulation parameter is always associated to a standard parameter
% A standard parameter can be associated to many modulation parameters
% I call extended parameter p_e a group consisting of a standard parameter
% and there associated modulation parameter
% p_e = p_s + ( p_m(1)*u(1)+...p_m(K)*u(K) )
% In the abscence of influence of u on p_s, p_m =0;
%
% Starting from a model as described above, and considering two possible
% states for the parametric modulation (actived/deactivated), the model
% defines and inverts all the combination of possible models.

% If there are N_m modulation parameters, there are 2^(N_m)
% models.
%
% Inversion is performed using data y and input u.
% A model comparison is performed and displayed
%-------------
% INPUT
% - y : vector of model output
% - u : vector of input to the model (experimenter data)
% - f_fname : handle of the evolution function
% - g_fname : handle of the observation function
% - dim : dimensions of the model
% - options : options
% - mod : structure declaring the parametric modulations
%   .(phi/theta).indp : matrix
%           Line 1 : index of the modulators (in param vector)
%           Line 2 : index of the modulating input (in input vector)
% - in : [irrelevant]
% OUTPUT
% - inversions : results of all inversions
%       .posterior
%       .out (see VBA_NLStateSpaceModel)
% - mod :
%   - .inv_order : binary matrix, each line describe wether or not modulation
% parameters are activated (1) or deactivated (0). modulation parameters
% are ordered as in output p_m
%   - .p_m : ordered cell of modulation parameters
%       .type : ('theta','phi') : class of parameter (evolution or observation)
%       .indp : index of the modulating parameter (in param vector)
%       .indu : index of the modulating input (in input vector)
% - LogEv : vector of the approximated log evidences of the ordered
% inversions.

posterior = [];
out = [];


N_p = 0;
try
for i_p = 1:length(mod.phi); % for each parameter
    N_p = N_p + length(mod.phi{i_p}.indpm);
end
for i_p = 1:length(mod.theta); % for each parameter
    N_p = N_p + length(mod.theta{i_p}.indpm);
end
catch
end

p_m = cell(1,N_p);


i_mod = 1;

try
for i_p = 1:length(mod.phi); % for each parameter
    for i_m = 1:length(mod.phi{i_p}.indu) % for each modulating input
        p_m{i_mod}.type = 'phi';
        p_m{i_mod}.indp = mod.phi{i_p}.indpm(1,i_m);
        p_m{i_mod}.indu = mod.phi{i_p}.indu(1,i_m);
        i_mod = i_mod+1;
    end
end
catch;
end

try
for i_p = 1:length(mod.theta); % for each parameter
    for i_m = 1:length(mod.theta{i_p}.indu) % for each modulating input
        p_m{i_mod}.type = 'theta';
        p_m{i_mod}.indp = mod.theta{i_p}.indpm(i_m);
        p_m{i_mod}.indu = mod.theta{i_p}.indu(i_m);
        i_mod = i_mod+1;
    end
end
catch;
end

mod.p_m = p_m;


% THis is where groups come into action!
try 
    for i_g = 1:length(mod.group_modulators)
        Ip = find(mod.group_modulators{i_g}.type == 'p');
        It = find(mod.group_modulators{i_g}.type == 't');
        mod.group_modulators{i_g}.indp_ordered = [Ip,It+dim.n_phi];
    end
catch;
    mod.group_modulators = cell(1,N_p);
    for i_mod = 1:N_p
        if  isequal(p_m{i_mod}.type, 'phi')
            mod.group_modulators{i_mod}.indp_ordered = i_mod;
            mod.group_modulators{i_mod}.type = 'p';
        elseif isequal(p_m{i_mod}.type, 'theta')
            mod.group_modulators{i_mod}.indp_ordered = i_mod + dim.n_phi;
            mod.group_modulators{i_mod}.type = 't';
        end
    end
end

% Ordering inversions
% All parameters are indexed : observation then evolution parameters

N_gm = length(mod.group_modulators); % number of groups of modulators

inv_order_gm = binary_cart_prod(N_gm);
N_inv = 2^N_gm;

if ~isequal(N_gm,N_p)
    inv_order = zeros(N_inv,N_p);
    for i_inv = 1:N_inv
        for i_gm = 1 : N_gm
            if inv_order_gm(i_inv,i_gm) == 1 % activate group
                inv_order(i_inv,mod.group_modulators{i_gm}.indp) = 1;
            end
            
        end
    end
else
    inv_order = inv_order_gm;
end

mod.inv_order = inv_order;
mod.inv_order_gm = inv_order_gm;


inversions = cell(1,N_inv); %all combinations
LogEv = zeros(1,N_p);

priors_default = options.priors;

for i_inv = 1:N_inv
    
    % load default priors
    priors = priors_default;
    
    % prepare priors for inversion
    for i_p = 1:N_p
        if (inv_order(i_inv,i_p) == 0)
            if isequal(p_m{i_p}.type,'theta')
                priors.muTheta(p_m{i_p}.indp) = 0;
                priors.SigmaTheta(p_m{i_p}.indp,p_m{i_p}.indp) = 0;
            elseif isequal(p_m{i_p}.type,'phi')
                priors.muPhi(p_m{i_p}.indp) = 0;
                priors.SigmaPhi(p_m{i_p}.indp,p_m{i_p}.indp) = 0;
            end
        end
    end
    
    options.priors = priors;
    [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
    
    inversions{i_inv}.out = out;
    inversions{i_inv}.posterior = posterior;
    
    
end


%-- Final display as proposed by Florent
LogEv = zeros(N_inv,1);
for i_inv = 1:N_inv
    LogEv(i_inv) =  inversions{i_inv}.out.F;
end


try modulator_names;
catch 
    modulator_names = cell(1,N_gm); 
end

mod.modulator_names = modulator_names;

VBA_Redisplay_modulation(LogEv,mod) % plot by groups


end

function o = binary_cart_prod(N_p)
% A binary count from 0 to 2^N_p-1 in the form a matrix
% Each line being the binary representation of a number
o = zeros(2^N_p,N_p);
for i = 0:2^N_p-1
    u = str2num(dec2bin(i)')';
    o(i+1,end-length(u)+1:end) =  u;
end
end


