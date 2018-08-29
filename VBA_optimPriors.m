function [optim_priors] = VBA_optimPriors(M,u,partition,density, priors2optim,Nsim)
% This functions performs the optimization over priors of models in order
% to minimize misclassifications within partitions of the set of models
% considered.
% E[e(y)] using either MCMC or Laplace
%
% INPUT
%  - M : structure containing descriptions of the structure of the models
%  considered (field options contains values of fixed priors)
%  - u : input to the models
%  - partition : cell of indices of models within each partition (models
%  are ordered)
%  - density : priors from which parameters are sampled to generate simulations and estimate the
%  predictive density
%  - Priors : structure containing for each model the grid of priors on
%  which optimization is performed
%  OUTPUT
%  - optim_priors : structure containing the optimal priors
%
%---------------------------------------------------------
% Definition of the structure of the variables used:
% - partition : cell 1*Nmodels
%       partition{p} : array 1*Nmodels_in_partition
% - M : cell 1*Nmodels
%       M{m} : structure
%           .f_fname
%           .g_fname
%           .options
%           .dim
% - density : cell 1*Nmodels
%       density{m} : structure (priors for simulations of model m)
% - Priors : cell 1*Nmodels
%       Priors{m}  : param2optim


%---------------------------------------------------------
%---  High level pseudo code
%---------------------------------------------------------
%
% Independent inversion of simulated data for each model
% - FOR each model that inverts the data
%   - FOR each configuration of the priors of this model
%     - FOR each model that generated the data
%       - FOR each simulation sample for this generative model
%           Compute the log evidence of the inversion
%         ENDFOR
%       ENDFOR
%     ENDFOR
%   ENDFOR
% Combination of all possible priors for all models
% Computation of the error
%
% - FOR each of these priors configurations
%   - FOR each generative model
%     - FOR each simulation from the generative model
%           Compute the error
%       ENDFOR
%       Sum errors over simulations
%     ENDFOR
%     Sum errors over generative models => this is the criteria for priors
%     optimization
%   ENDFOR
%
% Extracting the prior that maximized the error
%---------------------------------------------------------

%---------------------------------------------------------
% Defining the grid of priors on which to optimize & ordering simulations
% to be performed
%---------------------------------------------------------
% Here we restrict the definition of the grid by the definition of
% - the index of the prior to be optimized
% - the upper & lower bounds of the prior to be optimized
% - the grid size for each of the priors to be optimized


%----------------------------------------------------
% Check partitions
%----------------------------------------------------

disp('check partition.')
if isempty(partition) % create default partition (1 part = 1 model)
    partition = cell(1,numel(M));
    for p = 1 : numel(M)
        partition{p} = p;
    end
else % check if models are correctly indexed
    for p = 1 : numel(partition)
        for m = 1 : numel(partition{p})
            if ~all(ismember(partition{p},1:numel(M)))
                disp('partition is incorrect, optimization aborded')
                optim_priors = [];
                return
            end
        end
    end
end


% Running the simulations
disp('Running simulations from density...')
simulations = run_simulations(M,density,u,Nsim); % Simulations are made once for all
disp('done.')

% Ordering priors to be tested
disp('Preparing for optimization...')
priors_grid = prepare_optimization(M,priors2optim); % Organize grid exploration
disp('done.')

% Running the inversions of all data for each model
disp('Running inversions over grid...')
F = run_inversions(M,priors_grid,simulations,u,Nsim);
disp('done.')

% Computing the error for all combination of priors
disp('Computing error over grid...')
[e,V] = compute_error(M,F,priors_grid,Nsim,partition);
disp('done.')

% Extracting optimal priors
disp('Finding optima...')
optim_priors = return_optim_priors(M,e,V,priors_grid);
disp('done.')


optim_priors.partition = partition;

end


%-------------------------------------------------------------------------
% Functions
%-------------------------------------------------------------------------

%{
function e = getEfromY(y,u,M,m,k,gridp)
for m=1:length(M)
    opt = M{m}.options;
    opt.priors = Priors(priors,k,gridp)
    [posterior,out] = VBA_NLStateSpaceModel(y,u,M{m}.f_fname,M{m}.g_fname,M{m}.dim,opt)
    F(m) = out.F
end
e = scoreError(F,partition,m);
end
%}
%-------------------------------------------------------------------------

%{
function e = scoreError(F,partition,true)
F = F - max(F);
p = exp(F)./sum(exp(F));
for i=1:length(partition)
    pp(i) = sum(p(partition{i}));
end
e = 1 - pp(true);
end
%}
%-------------------------------------------------------------------------


function simulations = run_simulations(M,density,u,Nsim)
% This functions runs Nsim simulations of the models in M sampling
% parameters from their prior density
% INPUT
% - M : Structure of models
% - density : prior density on model parameters and initial states
% - u : input to models for simulation (shared for models and simulations)
% - Nsim : number of simulations per model
% OUTPUT
% - simulations : cell of size Nmodels containing output of simulations
%   .X : simulated hidden states (array)
%   .Y : simulated output (array)
Nmodels = length(M); % number of models considered
simulations = cell(1,Nmodels);
for m_gen = 1 : length(M)
    n_t = M{m_gen}.options.dim.n_t;  %length of time series
    opt = M{m_gen}.options; % simulation performed for model options
    opt.priors = density{m_gen}; % priors set to the chosen density for simulation
    
    try
    [pX,gX,pY,gY,X,Y,U] = VBA_MCMC_predictiveDensity_fb(M{m_gen}.f_fname,M{m_gen}.g_fname,u,n_t,opt,M{m_gen}.options.dim,Nsim,M{m_gen}.fb);
    simulations{m_gen}.X = X;
    simulations{m_gen}.Y = Y;    
    simulations{m_gen}.U = U;
    catch
    [pX,gX,pY,gY,X,Y] = VBA_MCMC_predictiveDensity(M{m_gen}.f_fname,M{m_gen}.g_fname,u,n_t,opt,M{m_gen}.options.dim,Nsim);
    simulations{m_gen}.X = X;
    simulations{m_gen}.Y = Y;
    end
end
end

%-------------------------------------------------------------------------

function priors_grid = prepare_optimization(M,priors2optim)
% This function orders the dimensions to be optimized for each model
% and orders the combinations of priors that will be evaluated
% This also performs a few checks
% INPUT
% - M : cell (1*Nmodels) containing the structure of the models
% - priors2optim : cell (1*Nodels) containing information about what priors
% should be optimized (bounds & steps on grid or grid values)
% OUTPUT
% - priors_grid : structure containing info about the ordering of priors
% combination for each models
%   - .ndim2opt (1*Nmodels)
%   - .ordered_priors2optim : cell of struct (1*Nmodels)
%       For each model, each dimension on which to optimize is given an index.
%       Indexed structure contains information about the prior:
%           - type : theta/phi
%           - moment : mu/s
%           - values : the array of values of the prior that will be
%           evaluated

Nmodels = length(M); % number of models considered
ndim2opt = zeros(1,Nmodels); % total number of dimension on which to optimize




%----------------------------------------------------
% Check optimization dimensions
%----------------------------------------------------

disp('----------Checking optimization dimensions------------')
for i_m = 1 : Nmodels % for each model
    disp(['---Model ',num2str(i_m)])
    dim_m = M{i_m}.options.dim;
    try
        fn = fieldnames(priors2optim{i_m}); %(x/theta/phi)
    catch
        fn = [];
    end
    for i = 1:numel(fn) % for each variable type
        fn1 = fieldnames(priors2optim{i_m}.(fn{i})); % (mu/s)
        for j = 1:numel(fn1) % for each moment
            disp([(fn{i}),'.',(fn1{j})])
            type = (fn{i});
            ind = priors2optim{i_m}.(fn{i}).(fn1{j}).ind;
            n = length(ind);
            ndim2opt(i_m) = ndim2opt(i_m)  + n;
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

%----------------------------------------------------
% Verify if priors given for priors that are not optimized are correct
%----------------------------------------------------

disp('----------Checking un-optimized priors ------------')

for i_m = 1:Nmodels
    
    disp(['---Model ',num2str(i_m)])
    priors = M{i_m}.options.priors;
    dim_m = M{i_m}.options.dim;
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
    M{i_m}.priors = priors;
end

%--------------------------------------
% Organizing dimensions
%--------------------------------------

disp('-------Organizing dimensions---------')

ordered_priors2optim = cell(1,Nmodels);
for i_m = 1:Nmodels
    ordered_priors2optim{i_m} = cell(1,ndim2opt(i_m));
end

for i_m = 1 : Nmodels % for each model
    i_dim = 0;
    disp(['---Model ',num2str(i_m)])
    
    try
        fn = fieldnames(priors2optim{i_m}); %(x/theta/phi)
    catch
        priors2optim{i_m} = [];
        fn = [];
    end
    for i = 1:numel(fn) % for each variable type
        fn1 = fieldnames(priors2optim{i_m}.(fn{i})); % (mu/s)
        
        for j = 1:numel(fn1) % for each moment
            
            disp([(fn{i}),'.',(fn1{j})])
            ind = priors2optim{i_m}.(fn{i}).(fn1{j}).ind;
            n = length(ind);
            
            for k = 1:length(ind)
                
                try % first try if vector of values was set
                    values = priors2optim{i_m}.(fn{i}).(fn1{j}).values(k,:);
                    
                catch % if it is not the case
                    
                    disp('Calculating grid values')
                    try % if it is not the case verify that bounds were set
                        
                        priors2optim{i_m}.(fn{i}).(fn1{j}).bounds(k,:);
                        
                        try % verify step
                            priors2optim{i_m}.(fn{i}).(fn1{j}).step(1,k);
                        catch % if not set it to 3
                            priors2optim{i_m}.(fn{i}).(fn1{j}).step(1,k) = 3;
                            disp('Warning : You have not set the grid size')
                            disp('set to default value: 3')
                        end
                        X = priors2optim{i_m}.(fn{i}).(fn1{j});
                        values = linspace(X.bounds(k,1),X.bounds(k,2),X.step(1,k));
                        
                    catch e
                        disp(e.message)
                        disp('Error : You have not set bounds for priors you wish to optim')
                        return;
                    end
                    
                    
                end
                
                i_dim = i_dim + 1;
                ordered_priors2optim{i_m}{i_dim}.model = i_m;
                ordered_priors2optim{i_m}{i_dim}.type = fn{i};
                ordered_priors2optim{i_m}{i_dim}.moment = fn1{j};
                ordered_priors2optim{i_m}{i_dim}.values = values;
                ordered_priors2optim{i_m}{i_dim}.ind = ind(k);
            end
        end
    end
end

priors_grid.ndim2opt = ndim2opt;
priors_grid.ordered_priors2optim = ordered_priors2optim;


%----------------------------------
% Computing all the combinations of parameters to be tested
%---------------------------------

disp('-------Preparing grid---------')
priors_combinations = cell(1,Nmodels);
Npriors_combinations = zeros(1,Nmodels);

for i_m = 1 : Nmodels % for each model
    if ~isempty(priors2optim{i_m})
        s = ['[',num2str(ordered_priors2optim{i_m}{1}.values),']'];
        for i = 2 : ndim2opt(i_m)
            s = [s,',[',num2str(ordered_priors2optim{i_m}{i}.values),']'];
        end
        eval(['priors_combinations{i_m} = cartprod(',s,');']);
        Npriors_combinations(1,i_m) = size(priors_combinations{i_m},1);
    else
        Npriors_combinations(1,i_m) = 1;
    end
end
priors_grid.Npriors_combinations = Npriors_combinations;
priors_grid.priors_combinations = priors_combinations;

end

%-------------------------------------------------------------------------

function F = run_inversions(M,priors_grid,simulations,u,Nsim)
% Performs the inversion of all simulations for all priors for each model.
% INPUT
% - M : cell of model structures and default priors
% - priors_grid : cell of structure containing info on priors optimization
% - simulations : cell of arrays of simulations of the different models
% - u : inputs used for simulations and inversion
% - Nsim : number of simulations per model
% OUTPUT
% - F : cell (1*Nmodels) of array of the results of inversions for all simulations

Nmodels = length(M);
F = cell(1,Nmodels);
Npriors_combinations = priors_grid.Npriors_combinations;

for m_inv = 1:length(M)
    
    disp(['Optimization for model ',num2str(m_inv),' out of ',num2str(Nmodels)])
    % for each 'inverting' model
    
    opt = M{m_inv}.options;
    opt.DisplayWin = 0;
    F_m = zeros(Npriors_combinations(m_inv),Nmodels,Nsim);
    
    for i_priors = 1:Npriors_combinations(m_inv)
        
        if mod(floor(i_priors/Npriors_combinations(m_inv)*100),10)==0
            disp([num2str(floor(i_priors/Npriors_combinations(m_inv)*100)),'%'])
        end
        
        % for each configuration of priors for the inverting model
        % opt.priors = getPriors(m_inv, i_priors,priors2optim);
        opt.priors = getPriors(opt.priors,m_inv,i_priors,priors_grid);
        for m_gen=1:length(M)
            %for each 'generating' model
            Y = simulations{m_gen}.Y;
            %load the N simulations for the generating model
            
            for i_sim=1:Nsim
                %for each of those simulations
                
                try
                    u = simulations{m_gen}.U(:,:,i_sim);
                catch
                end
                
                [posterior,out] = VBA_NLStateSpaceModel(Y(:,:,i_sim),...
                    u,...
                    M{m_inv}.f_fname,...
                    M{m_inv}.g_fname,...
                    M{m_inv}.options.dim,...
                    opt);
                
                %invert the model
                F_m(i_priors,m_gen,i_sim) = out.F;
                %and store the evidence
            end
            %e(i_sim) = scoreError(F(m_p,:,i_sim),partition,p);
            %from inversion of all models, I compute an error term
        end
        %Ee(p,m_p) = mean(e);
        % I approximate my expected error for each partition by a mean of errors over
        % samples
    end
    F{m_inv} = F_m;
    
end

end

%-------------------------------------------------------------------------

function priors = getPriors(priors,m_inv,i_priors,priors_grid)
% This loads the prior
% by modifying priors where optimized : priors
% of the model one wishes to invert : m_inv
% from the index : i_priors
% from prepared priors values : priors_grid
% INPUT
% - priors : standard structure, containing default priors
% - m_inv : index of model considered for inversion
% - i_priors : index of the prior configuration (as ordered in priors_grid)
% - priors_grid : cell of structs, ordered priors
% OUTPUT
% - priors : priors structure for the configuration i_priors
try
    p = priors_grid.priors_combinations{m_inv}(i_priors,:);
    for i_opt = 1 :length(p)
        type =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.type;
        moment =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.moment;
        ind =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.ind;
        if isequal(type,'phi')
            if isequal(moment,'mu')
                priors.muPhi(ind) = p(i_opt);
            elseif isequal(moment,'s')
                priors.SigmPhi(ind,ind) = p(i_opt);
            end
        elseif isequal(type,'theta')
            if isequal(moment,'mu')
                priors.muTheta(ind) = p(i_opt);
            elseif isequal(moment,'s')
                priors.SigmTheta(ind,ind) = p(i_opt);
            end
        end
    end
catch
    % unchanged priors
end
end

%-------------------------------------------------------------------------

function res = cartprod(varargin)
% Performs the cartesian product of arrays
% INPUT
% - varargin : arrays of values (1*N1), ... , (1*Nn)
% OUTPUT
% - res : array (N1*...*Nn)*Nn, cartesian product of input arrays
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

%-------------------------------------------------------------------------

function [E,index_of_combinations] = compute_error(M,F,priors_grid,Nsim,partition)
% Computes the mean of the probabilities that a model outside of the partition containing
% the generating models of simulations is identified as the model that
% actually generated the data.
% INPUT
% - M
% - F
% - priors_grid
% - Nsim
% OUTPUT
% - E
% - index_of_combinations

Nmodels = length(M);
Nparts = length(partition);
Npriors_combinations = priors_grid.Npriors_combinations;
s  ='';
for m = 1 : Nmodels
    s = [s,'[1:',num2str(max(Npriors_combinations(m),1)),']'];
    if m<Nmodels
        s = [s,','];
    end
end
eval(['index_of_combinations = cartprod(',s,');']);
Npriors_combinations_tot = size(index_of_combinations,1);

E = zeros(1,Npriors_combinations_tot);
for i_combi = 1 : Npriors_combinations_tot
    for p_gen = 1 : Nparts
        for m_gen = 1 : partition{p_gen}
            for i_sim = 1 : Nsim
                FF= [];
                for m_inv = 1 : Nmodels
                    FF(1,m_inv) = F{m_inv}( index_of_combinations(i_combi,m_gen),m_gen,i_sim);
                end
                FF = FF - max(FF);
                e(p_gen,i_sim) = sum(exp(FF(partition{p_gen}))/sum(exp(FF)));
                %e(m_gen,i_sim) = exp(FF(m_gen))/sum(exp(FF));
            end
        end
        E(i_combi) = sum(sum(e));
    end
end
end

%-------------------------------------------------------------------------

function optim_priors = return_optim_priors(M,e,V,priors_grid)
% Returns the optimal priors configurations in the form of a standard prior
% structure
% INPUT
% - M
% - e
% - V
% - priors_grid
% OUTPUT
% - optim_priors
Nmodels = length(M);
I = find(e == min(e)); % index of prior configuration that minimized the error
N_opt = length(I);
optim_priors.priors = cell(1,N_opt);

for i_min = 1:N_opt
    Priors = cell(1,Nmodels);
    % load values of optimized priors
    J = V(I(i_min),:); % indices of prior configuration per model
    for m_inv = 1 : Nmodels
        priors = M{m_inv}.options.priors;   % load default priors
        if priors_grid.priors_combinations{m_inv} > 0 % if necessary only
            p = priors_grid.priors_combinations{m_inv}(J(m_inv),:); % ordered priors values
            for i_opt = 1:length(p)
                type =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.type;
                moment =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.moment;
                ind =  priors_grid.ordered_priors2optim{m_inv}{i_opt}.ind;
                if isequal(type,'phi')
                    if isequal(moment,'mu')
                        priors.muPhi(ind) = p(i_opt);
                    elseif isequal(moment,'s')
                        priors.SigmPhi(ind,ind) = p(i_opt);
                    end
                elseif isequal(type,'theta')
                    if isequal(moment,'mu')
                        priors.muTheta(ind) = p(i_opt);
                    elseif isequal(moment,'s')
                        priors.SigmTheta(ind,ind) = p(i_opt);
                    end
                end
            end
        end
        Priors{m_inv} = priors;
    end
    optim_priors.priors{i_min} = Priors;
end

end


%-------------------------------------------------------------------------

