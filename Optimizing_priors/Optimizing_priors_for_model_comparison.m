function results = Optimizing_priors_for_model_comparison(simulations,models,optim_options)
% This script aims at optimizing priors of two models so that model
% selection based on evidence is most accurate.


%--------------------------------------
% IN
%   - simulations
%     1*Nmodels cell containing the simulations for each model
%         * 1*Nsim cell of struct containing each individual simulation
%           * x0
%           * theta
%           * phi
%           * u
%           * y
%           * isYout
%   - models
%     1*Nmodels cell
%         * 1*1 cell of struct containing model descriptions
%   - optim_options
%     - priors : as in classical inversion
%     - param2optim
%     1*Nmodels cell
%         *  1*1 cell of struct containing information on what to optimize
%     - method : not used yet : default is grid search
% OUT
%   - results
%--------------------------------------

%-------------
% The script, performs a basic grid search
% Inversions to be performed are ordered.
% From inversion index, one can find back what simulation was concerned and
% the parameters that were used.
%-------------

%------------- PSEUDO-CODE OF THE CLEF MODEL --------------%
% 1 - Verify optimisation dimensions
% 2 - Check unoptimized dimensions
%----------------------------------------------------------%

Nmodels = size(models,2); % Number of models that are compared
results = [];
%-----------------------------------------------------------------------
%----- STEP 1 : Choosing a criteria
%-----------------------------------------------------------------------

% In the case when there are only two models, we will for the moment chose
% a simple and single criteria :
% - For all simulated data, model selection must be as accurate as
% possible, that is the proportion of correctly selected models must be
% highest. This relies on a selection criteria that we take as having a difference of log evidence must be >3 in favor of the model that
% actually generated the data.
% Alternative criteria (more continuous ones) : one could be maximizing the
% oriented difference of log evidences

% INPUTS
% - criteria

%-----------------------------------------------------------------------
%----- STEP 2 : Optimizing the priors according the desired criteria
%-----------------------------------------------------------------------

% Many optimization schemes exsit.
% Among the possibilities, grid search or iterative algorithms.

% INPUTS
% - optimization scheme
% - stop criteria, maxIterations, etc
% - Subset of priors to consider for optimization
% - Range of these priors

% EXAMPLE : GRID SEARCH

%- Extracting the dimensions on which to optimize
%- counting dimensions on which to optimize





%--------------------------------------------------------------
% Verify if dimensions given for optimization match those of the models dimensions
%--------------------------------------------------------------

dim2opt = zeros(1,Nmodels); % total number of dimension on which to optimize

disp('----------Checking optimization dimensions------------')
for i_m = 1 : Nmodels % for each model
    
    disp(['---Model ',num2str(i_m)])
    dim_m = models{i_m}.options.dim;
    
    fn = fieldnames(optim_options.param2optim{i_m}); %(x/theta/phi)
    for i = 1:numel(fn) % for each variable type
        fn1 = fieldnames(optim_options.param2optim{i_m}.(fn{i})); % (mu/s)
        
        for j = 1:numel(fn1) % for each moment
            
            disp([(fn{i}),'.',(fn1{j})])
            
            type = (fn{i});
            ind = optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).ind;
            n = length(ind);
            dim2opt(i_m) = dim2opt(i_m)  + n;
            
            
            if isequal(type,'x')
                dim_m.n
                if ~isequal(ind <= dim_m.n,ones(size(ind)))
                    disp('optim dimensions exceed model dimension : hidden states')
                end
                
            elseif isequal(type,'theta')
                if ~isequal(ind <= dim_m.n_theta,ones(size(ind)))
                    disp('optim dimensions exceed model dimension : theta')
                end
                
            elseif isequal(type,'phi')
                if ~isequal(ind <= dim_m.n_phi,ones(size(ind)))
                    disp('optim dimensions exceed model dimension : phi')
                end
            else
            end
            
        end
    end
end


%--------------------------------------------------------------
% Verify if priors given for priors that are not optimized are correct
%--------------------------------------------------------------


disp('----------Checking un-optimized priors ------------')


for i_m = 1:Nmodels
    
    disp(['---Model ',num2str(i_m)])
    priors = optim_options.priors{i_m};
    dim_m = models{i_m}.options.dim;
    
    %------------------------------------------------
    try % ----------- muX0
        priors.muX0;
        if ~isequal(size(priors.muX0),[dim_m.n,1])
            disp('ERROR : incorrect dimension of muX0')
            return
        end
    catch
        disp('muX0 was not set -> set to 0')
        priors.muX0 = 0*ones(dim_m.n,1);
    end
    
    try  % ----------- SigmaX0
        priors.SigmaX0;
        if ~isequal(size(priors.SigmaX0),[dim_m.n,dim_m.n])
            disp('ERROR : incorrect dimension of SigmaX0')
            return
        end
    catch
        disp('SigmaX0 was not set -> set to 0')
        priors.SigmaX0 = 0*eye(dim_m.n);
    end
    %------------------------------------------------
    try % ----------- muTheta
        priors.muTheta;
        if ~isequal(size(priors.muTheta),[dim_m.n_theta,1])
            disp('ERROR : incorrect dimension of muTheta')
            return
        end
    catch
        disp('muTheta was not set -> set to 0')
        priors.muTheta = 0*ones(dim_m.n_theta,1);
    end
    
    try  % ----------- SigmaTheta
        priors.SigmaTheta;
        if ~isequal(size(priors.SigmaTheta),[dim_m.n_theta,dim_m.n_theta])
            disp('ERROR : incorrect dimension of SigmaTheta')
            return
        end
    catch
        disp('SigmaTheta was not set -> set to 0')
        priors.SigmaTheta = 0*eye(dim_m.n_theta);
    end
    %------------------------------------------------
    try % ----------- muPhi
        priors.muPhi;
        if ~isequal(size(priors.muPhi),[dim_m.n_phi,1])
            disp('ERROR : incorrect dimension of muPhi')
            return
        end
    catch
        disp('muPhi was not set -> set to 0')
        priors.muPhi = 0*ones(dim_m.n_phi,1);
    end
    
    try  % ----------- SigmaPhi
        priors.SigmaPhi;
        if ~isequal(size(priors.SigmaPhi),[dim_m.n_phi,dim_m.n_phi])
            disp('ERROR : incorrect dimension of SigmaPhi')
            return
        end
    catch
        disp('muPhi was not set -> set to 0')
        priors.SigmaPhi = 0*eye(dim_m.n_phi);
    end
    
    optim_options.priors{i_m} =priors;
    
    
end

%--------------------------------------------------------------
% Organizing dimensions

disp('-------Organizing dimensions---------')

dim_optim = cell(1,Nmodels);
for i_m = 1:Nmodels
    dim_optim{i_m} = cell(1,dim2opt(i_m));
end

for i_m = 1 : Nmodels % for each model
    i_dim = 0;
    disp(['---Model ',num2str(i_m)])
    
    
    fn = fieldnames(optim_options.param2optim{i_m}); %(x/theta/phi)
    for i = 1:numel(fn) % for each variable type
        fn1 = fieldnames(optim_options.param2optim{i_m}.(fn{i})); % (mu/s)
        
        for j = 1:numel(fn1) % for each moment
            
            disp([(fn{i}),'.',(fn1{j})])
            ind = optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).ind;
            n = length(ind);
            
            for k = 1:length(ind)
                
                try % first try if vector of values was set
                    values = optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).values(k,:);
                    
                catch % if it is not the case
                    
                    disp('Calculating grid values')
                    try % if it is not the case verify that bounds were set
                        
                        optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).bounds(k,:);
                        
                        try % verify step
                            optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).step(1,k);
                        catch % if not set it to 3
                            optim_options.param2optim{i_m}.(fn{i}).(fn1{j}).step(1,k) = 3;
                            disp('Warning : You have not set the grid size')
                            disp('set to default value: 3')
                        end
                        X = optim_options.param2optim{i_m}.(fn{i}).(fn1{j});
                        values = linspace(X.bounds(k,1),X.bounds(k,2),X.step(1,k));
                        
                    catch e
                        disp(e.message)
                        disp('Error : You have not set bounds for priors you wish to optim')
                        return;
                    end
                    
                    
                end
                
                i_dim = i_dim + 1;
                dim_optim{i_m}{i_dim}.model = i_m;
                dim_optim{i_m}{i_dim}.type = fn{i};
                dim_optim{i_m}{i_dim}.moment = fn1{j};
                dim_optim{i_m}{i_dim}.values = values;
                dim_optim{i_m}{i_dim}.ind = ind(k);
                
            end
            
            
            
        end
    end
end

results.dim2opt = dim2opt;
results.dim_optim = dim_optim;

%-------------------------------------------------------------------------------
% Computing all the combinaisons of parameters to be tested
%-------------------------------------------------------------------------------
disp('-------Preparing grid---------')

combi = cell(1,Nmodels);


for i_m = 1 : Nmodels % for each model
    
    s = ['[',num2str(dim_optim{i_m}{1}.values),']'];
    for i = 2 : dim2opt(i_m)
        s = [s,',[',num2str(dim_optim{i_m}{i}.values),']'];
    end
    
    
    eval(['V = cartprod(',s,');']);
    combi{i_m} = V;
    
    results.combi = combi;
    
end

%-------------------------------------------------------------------------------
% Performing the inversions
%-------------------------------------------------------------------------------
disp('-------Performing inversions---------')




%--------- preparing order of simulations

Nsim = 0;
Nsim_per_model = [];
sim = []; % index of generating model, index of sim per generative model
for i_m = 1 : Nmodels
    n =size(simulations{i_m},2);
    Nsim_per_model(i_m) = n;
    sim = [sim;[i_m*ones(n,1),(1:n)']];
end
Nsim = sum(Nsim_per_model);



Ninv_per_model=[];
inv = cell(1,Nmodels);
for i_m = 1 : Nmodels % for all models
    Ninv_per_model(1,i_m) = size(combi{i_m},1);
    inv{i_m} = cell(1,Ninv_per_model(1,i_m));
end


% inversions are performed separately
% - for each model
%    - for each simulation
%       - for each possible set of priors (cartesian product)
%         end
%      end
%   end




for i_m_inv = 1 : Nmodels % for all models
    
    c = combi{i_m_inv};
    Ninv_per_model = size(c,1);
    
    
    f_fname = models{i_m_inv}.f_fname;
    g_fname = models{i_m_inv}.g_fname;
    dim = models{i_m_inv}.options.dim;
    
    %--------- Setting parameters that are not optimized
    %--------- (common to all inversions)
    
    options = models{i_m_inv}.options;
    priors = optim_options.priors{i_m_inv};
    
    for i_m_sim = 1 : Nmodels
        
        for i_sim = 1 : Nsim_per_model(i_m_sim) % for all simulations
            i_sim = sim(i_sim,2); % index of sim within model category
            
            %--------- loading simulation data
            i_sim
            
            u = simulations{i_m_inv}{i_sim}.u;
            y = simulations{i_m_inv}{i_sim}.y;
            options.isYout = simulations{i_m_inv}{i_sim}.isYout;
            
            
            for i_inv = 1 : Ninv_per_model(1,i_m_inv) % for all predefined priors
                
                disp(['Inverted sim #',num2str(i_inv)])
                
                %--------- Setting prior parameters that have to be optimized
                for i_p = 1 : dim2opt(i_m_inv)
                    
                    ind = dim_optim{i_m_inv}{i_p}.ind;
                    type = dim_optim{i_m_inv}{i_p}.type;
                    moment = dim_optim{i_m_inv}{i_p}.moment;
                    values = dim_optim{i_m_inv}{i_p}.values;
                    
                    if isequal(type,'x0')
                        if isequal(moment,'mu')
                            priors.muX0(1,ind) = c(i_inv,i_p);
                        elseif isequal(moment,'s')
                            priors.SigmaX0(ind,ind) = c(i_inv,i_p); % diag term
                        end
                        
                    elseif isequal(type,'theta')
                        if isequal(moment,'mu')
                            
                            priors.muTheta(1,ind) = c(i_inv,i_p);
                        elseif isequal(moment,'s')
                            priors.SigmaTheta(ind,ind) = c(i_inv,i_p); % diag term
                        end
                        
                    elseif isequal(type,'phi')
                        if isequal(moment,'mu')
                            priors.muPhi(1,ind) = c(i_inv,i_p);
                        elseif isequal(moment,'s')
                            priors.SigmaPhi(ind,ind) = c(i_inv,i_p); % diag term
                        end
                    end
                    
                end
                
                options.priors = priors;
                
                %--------- Performing inversion ------------
                
                disp('oups')
                %  [posterior,out] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
                
                %--------- Saving inversion's result ------------
                %  inv{i_m}{i_inv}.F = out.F;
                % inv{i_m}{i_inv}.inv = i_m;
                % inv{i_m}{i_inv}.sim = sim(i_sim,1);
                
                
            end
            
        end
    end
    
end


%-------------------------------------------------------------------------------
% Finding the optimal priors
%-------------------------------------------------------------------------------

% until now, inversions were performed separately for each model.
%--- Criteria here :  sum of log evidences

% FOCUS on the cas of 2 models
% Here we want to
% - maximize the evidence of a model for simulations it generated
% - minimize the evidence of a model for simulations it didn't
% Doing so will maximize the difference of evidence in favor of the model
% that actually generated the data


% Here is where the crossing takes place

% the thing should go that way
% - for each set of parameters
%  compute the sum on


% for i_m = 1 : Nmodels
%
%     i_max = 0;
%     F_max = -Inf;
%
%     for i_inv = 1 : Ninv(i_m) %
%
%         F=inv{i_m}{i_inv}.F;
%         sim = inv{i_m}{i_inv}.sim;
%
%
%         if (sim == i_m) % Model that INVERTS = Model that GENERATES
%             if F>F_max
%                 F_max = F;
%                 i_max = i_inv;
%             end
%
%         else % Model that INVERTS =/= Model that GENERATES
%             if F>F_max
%                 F_max = F;
%                 i_max = i_inv;
%             end
%
%         end
%
%
%     end
%
% end




end


function res = cartprod(varargin) % cartesian product
res = [];
if nargin >= 1
    x = varargin{1};
    x = x(:);
    if nargin == 1
        res = x;
        return;
    end
    y = varargin(2:end);
    y = cartprod(y{:});
    
    nx = size(x, 1);
    ny = size(y, 1);
    rx = repmat(x.', ny, 1);
    ry = repmat(y, nx, 1);
    res = [rx(:), ry];
end
end
