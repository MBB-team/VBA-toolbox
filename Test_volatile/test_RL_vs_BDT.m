% In this script, I simulate a to models of learning and decision making on
% an operant conditioning task

clear all
close all
clc

%----------------------------------------------------------------------
% Description of the task:
% - 2 alternatives
% - for each alternative, 2 possible outcomes with different probabilities
% What can be manipulated in the task
% - magnitude of each outcome
% - probability of each outcome for each alternative


%----------------------------------------------------------------------
% Different conditions of interest :
%----- Probabilities :
% 4 conditions : 25/75 ; 75/25 ; 75/75 ; 25/25
%----- Values :
% 3 conditions : +/0 ;  -/0 ; -/+
% Condition 1 : o1 = 10; o2 = 0;
% Condition 2 : o1 = -10; o2 = 0;
% Condition 3 : o2 = -10; o2 = 10;
%----------------------------------------------------------------------



OUTCOMES = [10 , 0];
    %...
    %-10, 0;
    %-10, 10];

PROBABILIES = [.40, .60];
    %...
    %.75, .75; ...
    %.75, .25;
    %.25, .25];

No = size(OUTCOMES,1);
Np = size(PROBABILIES,1);

Ntrials = 200;
Nit = 1;

figure(3)

for i_o = 1:size(OUTCOMES,1)
    for i_p = 1:size(PROBABILIES,1);
                
        % in order to plot the learning curve, we need to save decisions
        DecisionsBDT = zeros(Nit,Ntrials);
        LearnedProba = zeros(Nit,Ntrials,2);
        
        
        for it = 1 : Nit
               
        p1_1 = PROBABILIES(i_p,1); % probability of outcome 1 for choice 1
        p1_2 = PROBABILIES(i_p,2); % probability of outcome 1 for choice 2
        o1 = OUTCOMES(i_o,1); % value of outcome 1
        o2 = OUTCOMES(i_o,2); % value of outcome 2            
        
        condtxt = ['proba:',num2str(p1_1),'/',num2str(p1_2),', outcomes:',num2str(o1),'/',num2str(o2)];
        
%         p1_1 = 0.9; % probability of outcome 1 for choice 1
%         p1_2 = 0.2; % probability of outcome 1 for choice 2
%         o1 = 10; % value of outcome 1
%         o2 = 2;  % value of outcome 2
        
        % Generating rewards for each alternatives
        i1 = rand(1,Ntrials)<p1_1; % outcome for choice 1
        i2 = rand(1,Ntrials)<p1_2; % outcome for choice 2
        R = [o1.^i1.*o2.^(1-i1);o1.^i2.*o2.^(1-i2)];
        
        %----------------------------------------------------------------------
        % Description of the competing models
        
        %%
    %{
    %----------------- RL model + softmax decision rule
        
        % For any predefined rewards for all alternatives
        h_fname = @h_reward_2Q;
        fb.inH.u0 = R; % with reversals
        fb.h_fname = h_fname;
        fb.indy = 1;
        fb.indfb = [2,3]; % indices where put feedbacks in the experimenter data matrix u
        
        %%% Comparison to what I had done
        % feedbacks where specified in advance in u. And simulation just added the
        % the generated data in the first line.
        % here feedbacks are specified in the feedback structure
        
        alpha = 0.2; % learning rate
        beta =0.5; % inverse temperature
        Q0 = [0;0];
        
        f_fname = @f_Qlearn_2Q;
        g_fname = @g_softmax_2Q;
        theta = sigm(alpha,struct('INV',1));
        phi = log(beta);
        u = [zeros(1,Ntrials);R];
        alpha = Inf;
        sigma =  Inf;
        
        dim_output = 1; % the delta reaction-time
        dim_data = 3; % index of sequence(not used for sim)/ rewards for both alternatives
        dim = struct('n',2,...  %( 2 (Qvalues) * Nsessions
            'p',dim_output,... % total output dimension
            'n_theta',1,... % evolution parameters
            'n_phi', 1,... % observation parameters
            'n_t',Ntrials);
        
        options.inF = [];
        options.inG = [];
        options.DisplayWin = 1;
        options.GnFigs = 0;
        options.binomial = 1; % Dealing with binary data
        options.isYout = zeros(1,Ntrials); % Excluding data points
        options.dim = dim;
        
        options.skipf = zeros(1,Ntrials);
        options.skipf(1) = 1; % apply identity mapping from x0 to x1.
        
        [y,x,x0,eta,e,u] = simulateNLSS_fb(Ntrials,f_fname,g_fname,theta,phi,u,Inf,Inf,options,Q0,fb);
        
        DecisionsRL(it,:) = y; % saving choice sequence
        
        %figure(1)
        %plot(x(1,:))
        %title('Qvalues of each alternatives')
   %}     
        %%
        %----------------- Probability learning + utility + softmax
        
        
        
        % Remarks on structure :
        % - 2 models based on Mathys, Daunizeau 2010, tracking outcome probability
        % - 5 hidden states each
        
        
        % evolution, observation and feedback functions
        f_fname = @f_OpLearn; % dual Mathys with feedback
        %g_fname = @g_VBvolatile0; % evolution function : softmax + bias (Jean)
        %g_fname = @g_VBvolatile1; % Probability matching
        g_fname = @g_VBvolatile2; % utility + Softmax
        h_fname = @h_truefalse; % returns reward (0 or 1)
        %h_fname = @h_rewardpunish; % returns reward (-1 or 1)
        
        % defining the utility function
        u_fname = @u_prospect_theory; % handle of utility function
        inG.u_fname = u_fname;
        inG.o1 = o1;
        inG.o2 = o2;
        % allocate feedback struture for simulations
        
        fb.inH.u0 = R; % with reversals
        fb.h_fname = @h_reward_2Q;
        fb.indy = 1;
        fb.indfb = [2,3];
        
        
        % simulation parameters
        inF.lev2 = 1; % remove 3rd level (volatility learning)
        inF.kaub = 1.4;
        inF.thub = 1;
        inF.rf = -1;
        inG.respmod = 'taylor';
        
        % choose initial conditions
        x0 = repmat([0.5;0;0;1;log(4)],2,1); % identical initial conditions for both models
        %u = zeros(2,size(fb.inH.u0,2)+1);
        u = [zeros(1,Ntrials);R];
        
        dim_output = 1; % the delta reaction-time
        dim = struct('n',2*5,...
            'p',dim_output,... % total output dimension
            'n_theta',3,...
            'n_phi',2,...
            'n_t',Ntrials);
        
        % phi: parameter of utility function
        % theta : kappa,omega,theta
        
        theta = [1;-4;-1];
        phi = [log(1);log(1)]; % inverse temperature in softmax, parameter of utility function
        
        priors.muPhi = phi;% zeros(dim.n_phi,1);
        priors.muTheta = [0;-4;0];
        priors.muX0 = x0;
        priors.SigmaPhi = 1e2*eye(dim.n_phi);
        priors.SigmaTheta = 1e2*eye(dim.n_theta);
        priors.SigmaX0 = 0e1*eye(dim.n);
        priors.a_alpha = Inf;
        priors.b_alpha = 0;
        
        options.priors = priors;
        options.binomial = 1;
        options.inF = inF;
        options.inG = inG;
        options.dim = dim;
        
        options.skipf = zeros(1,length(u));
        options.skipf(1) = 1; % apply identity mapping from x0 to x1.
        
        [y,x,x0,eta,e,u] = simulateNLSS_fb(length(u),f_fname,g_fname,theta,phi,u,Inf,Inf,options,x0,fb);
        DecisionsBDT(it,:) = y; % saving choice sequence
        LearnedProba(it,:,1) = sigm(x(2,:));
        LearnedProba(it,:,2) = sigm(x(7,:));

        %{
        figure(2)
        subplot(1,2,1)
        plot(y-e,'r')
        hold on
        plot(y,'kx')
        legend({'p(y=1|theta,phi,m)','binomial data samples'})
        getSubplots
        
        subplot(1,2,2)
        hold on
        plot(x(1,:),'r')
        plot(x(2,:),'--r')
        
        plot(x(6,:),'b')
        plot(x(7,:),'--b')
        %}
        
        
       
        end
        
        
        
        %-------------- Plotting learning curves
        % Proportion of choice of the most rewarding task.
        
        figure(4)  % BDT : Learned probabilities
        hold on     
        subplot(No,Np,(i_o-1)*Np+i_p)
        hold on
        plot(mean(LearnedProba(:,:,1),1))
        plot(mean(LearnedProba(:,:,2),1))
        title([condtxt])

        
        
        figure(3)  % Bayesian Decision Theory : LEARNING CURVES
        hold on     
        subplot(No,Np,(i_o-1)*Np+i_p)
        lc = @(a,x)(0.5 + a(1)*(1-exp(-a(2)*x)));
        %beta = nlinfit(1:Ntrials,mean(DecisionsBDT,1),lc,[0,1]);
        hold on
        plot(mean(DecisionsBDT,1))
        %plot(1:Ntrials,lc(beta,1:Ntrials),'r');
        title([condtxt])

    
       
        
        
        
        
        
    end
end

%%
% figure
% hold on
% plot(posterior.muX(2,:),'r')
% plot(posterior.muX(2+5,:),'b')
% plot(y,'x')
%
% figure
% hold on
% plot(posterior.muX(1,:),'r')
% plot(posterior.muX(6,:),'b')
% plot(y,'x')
%
% title('bern')
