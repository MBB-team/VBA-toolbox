function [posterior,out] = demo_multisource()
% This script demonstrate how to simulate and invert multisource models, 
% ie using multiple data following different distributions
%
% In this demo, we use a model which predict simultaneously a two
% dimensional gaussian observation and a binomial sequence

%% === model specification =======================================
% the observation function predict a vector of value corresponding to the 
% vertical concatenation of all sources (gaussians and binomial)

    function [gx] = g_demo_multisource(Xt,P,ut,in) 
        % linear mapping for the normal data
        gx(1) = P(1) + P(2)*ut(1) ;   
        gx(2) = P(1) + P(3)*ut(2) ;   
        % logisitic mapping for the binomial data
        gx(3) = 1./(1+exp(- P(4) - P(5)*(ut(1)+ut(2)))) ;   
        % insure verticality
        gx = vec(gx);
    end

%% === data distribution =======================================
% the field 'sources' of the 'options' structure explicits how to split the
% observations into the different distribution types

sources(1) = struct('out',1,'type',0); % the first two lines of data are gaussian with independent variance 
sources(2) = struct('out',2,'type',0); % set sources(1) = struct('out',1:2,'type',0) to use common hyperparameters
sources(3) = struct('out',3,'type',1); % the third line of data is binomial 

options.sources=sources;

%% === simulate data =======================================

% generate experimental setup
T = 200;
u = randn(2,T);

% Choose basic settings for simulations
g_fname = @g_demo_multisource; % observation function
phi = [1;10;15;-2;4];          % observation parameters
sigma = [1e-1 1e0];            % respective precision of the gaussian sources


disp('*** Simulation');
[y,x,x0,eta,e] = simulateNLSS(T,[],g_fname,[],phi,u,Inf,sigma,options,[]);


%% === invert model =======================================
        
dim.n_phi = 5;                  % nb of observation parameters
dim.n_theta = 0;                % nb of evolution parameters
dim.n=0;                        % nb of hidden states

[posterior,out] = VBA_NLStateSpaceModel(y,u,[],g_fname,dim,options);

end