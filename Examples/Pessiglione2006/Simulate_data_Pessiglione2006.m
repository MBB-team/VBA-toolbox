function [Y,U,IsYout] = Simulate_data_Pessiglione2006()




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Description of the experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsubject = 3;
% 3 groups
N1 = 1; % placebo (single blind)
N2 = 1; % (Haldol)
N3 = 1; % (levodopa)
group = [ones(1,N1),ones(1,N2)*2,ones(1,N3)*3];

% 3 independant sessions per subject
Nsessions_ps = 3;
% dividing each session in 3 (for each condition) 
Npairs = 2;
% leads to 3*3=9 independent sessions per subject
Nsessions = Npairs*Nsessions_ps;

Nsessions_tot=Nsessions*Nsubject;

NtrialsPerSession = 15*ones(1,Nsessions_tot); % 30
Ntrials = max(NtrialsPerSession); 


% Generating parameters according to reported fits
% The constants ? (learning rate) and ? (temperature) were adjusted to maximise 
% the probability (or likelihood) of the actual choices under the model. 
% For the gain / loss conditions respectively, we found 
% ? = 0.29 / 0.46 and ? = 0.18 / 0.33, 
% with 95% confidence intervals of 0.24-0.31 / 0.40-0.52 and 0.17-0.20 / 0.31-0.35

alpha_gain = 0.29; % IC95% 0.24-0.31 
alpha_loss = 0.46; % IC95% 0.40-0.52
beta_gain = 0.18; % IC95% 0.17-0.20
beta_loss = 0.33; % IC95% 0.31-0.35
alpha = 0.2;%[0.29,0.46,0];
beta = 5;%1./[0.18,0.33,0.2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Simulating behavioral data and formating data for inversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% For each pair of cue we order them so that they give OUTCOME or 0 with
% probability PROB
OUTCOME = [1   ,-1  ,0  ];
PROB    = [0.8 ,0.8 ,0.8;
           0.2 ,0.2 ,0.2];


dim_output = 1; % number of output data per session (only choice response)
dim_data = 4; % number of experimenter data per session  
%1:action of subject
%2:reward obtained by subject
%3:which group of subject
%4:pair of cues (gain/loss)


U = zeros(Nsubject*Nsessions_ps*Npairs*dim_data,Ntrials);
Y = zeros(Nsubject*Nsessions_ps*Npairs,Ntrials);
IsYout = Y;



for subj = 1 : Nsubject % for each subject

    for s = 1 : Nsessions_ps % for each session
        
        for p = 1 : Npairs % for each of the 3 pairs of cues
        
         % Generate rewards for each alternative   
         R = (rand(2,Ntrials)<PROB(:,p)*ones(1,Ntrials))*OUTCOME(p);
         % Initial Qvalues
         Q0 = [0;0];
         
         % Simulate QLearning + softmax decisions
         [A,r,Q] = simulate_QLearning_2Q(Q0,alpha(1),beta(1),R);
         
         % Concatenate data and output for further inversion
         
         i = Nsessions*(subj-1) + Npairs*(s-1) + p;
         
         U((i-1)*dim_data+1:(i)*dim_data,1:NtrialsPerSession(i)) = [A;r;group(subj)*ones(1,Ntrials);p*ones(1,Ntrials)];
         Y((i-1)*dim_output+1:(i)*dim_output,1:NtrialsPerSession(i)) = A;
         IsYout((i-1)*dim_output+1:(i)*dim_output,1:NtrialsPerSession(i)) = 0;
         
         
         
        end
        
    end 
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% some plots of performance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Separating conditions

close all

AG = U(1:9:Nsessions_tot,:);
AL = U(4:9:Nsessions_tot,:);

T = 1:Ntrials;
figure
hold on
errorfill(T,1-mean(AG),std(AG),'r')
errorfill(T,mean(AL),std(AL),'g')
title('proba of choice of the Correct option (red:GAIN/green:LOSS)')
xlabel('Trials')
