function [inversions,inv_order,p_m,LogEv] = VBA_Inversion_modulation(y,u,f_fname,g_fname,dim,options,mod,in)
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
% - inv_order : binary matrix, each line describe wether or not modulation
% parameters are activated (1) or deactivated (0). modulation parameters
% are ordered as in output p_m
% - p_m : ordered cell of modulation parameters
%       .type : ('theta','phi') : class of parameter (evolution or observation)
%       .indp : index of the modulating parameter (in param vector)
%       .indu : index of the modulating input (in input vector)
% - LogEv : vector of the approximated log evidences of the ordered
% inversions.

posterior = [];
out = [];


% Ordering

try mod.phi.indp;
catch
    mod.phi.indp = []; end

try  mod.theta.indp;
catch
    mod.theta.indp = []; end


% Ordering modulation parameters
N_p = size(mod.phi.indp,2)+size(mod.theta.indp,2);
p_m = cell(1,N_p);
modulator_names = cell(N_p,1);
i_p = 0;
for indp = mod.theta.indp
    i_p = i_p+1;
    p_m{i_p}.type = 'theta';
    p_m{i_p}.indp = indp(1);
    p_m{i_p}.indu = indp(2);
    
    imp =  indp(1);% index of modulated parameter
    imu = mod.indu(indp(2));% index of modulating input
    modulator_names{i_p} = ['u(',num2str(imu),')->theta(',num2str(imp),')'];
end

for indp = mod.phi.indp
    i_p = i_p+1;
    p_m{i_p}.type = 'phi';
    p_m{i_p}.indp = indp(1);
    p_m{i_p}.indu = indp(2);
    
    imp =  indp(1);% index of modulated parameter
    imu = mod.indu(indp(2));% index of modulating input
    modulator_names{i_p} = ['u(',num2str(imu),')->phi(',num2str(imp),')'];
    
end

% Ordering inversions

% THis is where groups come into action! after poo
try mod.group_modulators;
catch;
    mod.group_modulators = cell(1,N_p);
    for i = 1:N_p
        mod.group_modulators{i}.indp = i;
    end
end
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


%plot_mod(LogEv,inv_order,modulator_names)
VBA_Redisplay_modulation(LogEv,inv_order,modulator_names)

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

function [] = plot_mod(LogEv,inv_order,modulator_names)

Nsubject = size(LogEv,2); % number of subjects considered

if Nsubject == 1 % Case single subject
    N_p = size(inv_order,1);
    
    Nmodels = N_p;
    
    pp = exp(LogEv-max(LogEv));
    pp = pp./sum(pp);
    
    figure
    
    subplot(1,3,1)
    imagesc(-inv_order)
    colormap(gray)
    
    if ~isempty(modulator_names)
        
        N_m = size(modulator_names,1); % number of modulators
        set(gca, 'XTickLabel',modulator_names) % setting label names
        set(gca,'XLim',[0.5 N_m+0.5])
        set(gca,'XTick',[1:N_m]) % setting label positions
        
        set(gca,'XTickLabel',[]);%erase current tick labels from figure
        c=get(gca,'YTick')+N_p-1;%make new tick labels
        b=get(gca,'XTick');%get tick label positions
        text(b,repmat(N_p+1,N_m,1),modulator_names,'HorizontalAlignment','right','rotation',45); % rotating labels
        
    else
        xlabel('Modulation parameters', ...
            'FontSize', 10, ...
            'FontWeight', 'bold')
    end
    
    
    ylabel('Models')
    
    subplot(1,3,2)
    barh((LogEv-max(LogEv)),0.5)
    axis([min((LogEv-max(LogEv))) max((LogEv-max(LogEv))) 1-0.5 Nmodels+0.5])
    set(gca, 'YDir', 'reverse')
    xlabel('Log Evidence')
    
    subplot(1,3,3)
    barh(pp,0.5)
    axis([0 max(pp) 1-0.5 Nmodels+0.5])
    set(gca, 'YDir', 'reverse')
    xlabel('Posterior probability')
    
    
else % Case multiple subject
    
    N_p = size(inv_order,1);
    Nmodels = 2^N_p;
    
    pp = exp(LogEv-max(LogEv));
    pp = pp./sum(pp);
    
    figure
    subplot(1,3,1)
    imagesc(-inv_order)
    colormap(gray)
    xlabel('Modulation parameters')
    ylabel('Models')
    
    
    S_LogEv = sum(LogEv,2); % sum of log-evidence
    subplot(1,3,2)
    barh((S_LogEv-max(S_LogEv)),0.5)
    axis([min((S_LogEv-max(S_LogEv))) max((S_LogEv-max(S_LogEv))) 1-0.5 Nmodels+0.5])
    set(gca, 'YDir', 'reverse')
    xlabel('Log Evidence')
    
    % performing fixed effects analysis
    [exp_r,xp,r_samp,g_post] = spm_BMS_gibbs (LogEv, ones(1,Nmodels))
    
    subplot(1,3,3)
    barh(xp,0.5)
    axis([0 1 1-0.5 Nmodels+0.5])
    set(gca, 'YDir', 'reverse')
    xlabel('Exceedance probability')
    
    
    
end

end

