function data_tache_mere = LOAD_and_FORMAT_data_tache_mere(  )

% This scripts format the data of the tache mère in order to easily perform
% inversion

clear all
close all


%-----NOTES on how file is organized
% session1 = training session
% session2-3 = test sessions
%-- DATA
% column 1 = index of session
% column 2 = index of trial
% column 3 = condition (1=Gain (1€), 2=Neutral (0€), 3=Loss (-1€))
% column 4 = checktime
% column 5 = 1=Go, -1=NoGo
% column 6 = choice (1=correct, -1=incorrect)
% (correct = most rewarded option in cond=1 | less punished option in cond=3)
% colonne 7 = feedback (cond 1 : 1 = 1€ , -1=0€ | cond 3, 1=0€, -1=-1€)
% colonne 8 = reaction time, (millisec if Go, 0 if NoGo)
% colonne 9 = number of hesitations



%----------- DO THE INVERSION FOR ALL SUBJECTS

controls1=[1001:1015 1103 1106 1108 1109]; % controles agés (Stefano, tache mêre)
data_tache_mere = cell(1,length(controls1)); % where results for each subjects will be saved
SESSIONS = [2,3]; % the two training sessions

i_s =0; % index of subject
for i_subject = controls1
    
    i_s = i_s+1;
    %----------- FIRST : Loading the data & formating for inversion
    
    root = 'C:\Users\vincent.adam\Desktop\Vincent ADAM\Matlab\Demo_toolbox\Data\Data_tache_mere\VINCENT\'; % location of data folder
    subfolder = 'CONTROLS_1\';
    
    %-- LOADING THE TWO SESSIONS --%
    
    Y_PL = []; Y_QL = [];
    U_PL = []; U_QL = [];
  
    Outcomes = [1,0;... % Gain condition
                0,-1];  % Loss condition
    for i_session = SESSIONS; % index of session
        
        filename = ['ILtestSub',num2str(i_subject),'Session',num2str(i_session)]; % as saved by stefano
        file_root = [root, subfolder, filename]; % position of file on hard-drive
        
        load(file_root); % loading file        
        
        %---- Formating data for further inversion
        
        cond = data(:,3); % vector of conditions
        i_cond1 = find(cond==1); % indices of trial where condition is 1 (gain)
        i_cond3 = find(cond==3); % indices of trial where condition is 3 (loss)
        
        y1 = data(i_cond1,6)'>0; % binary choice (0: correct, 1:incorrect)
        y3 = data(i_cond3,6)'>0; % binary choice (0: correct, 1:incorrect)
        
        %--------------------------------------------------------------
        % FORMATING DATA FOR PROBABILITY LEARNING
        %--------------------------------------------------------------
        % NOTES on how data should be formated
        % y: binary choice
        % u:
        %   - subject choice (binary)
        %   - outcome category (binary)
        
        % outputs
          % experimenter data
        u1_PL = [y1;... % binary choice (0: correct, 1:incorrect)
            data(i_cond1,7)'>0]; % binary outcome (0: best, 1:worst)
        u3_PL = [y3;... % binary choice (0: correct, 1:incorrect)
            data(i_cond3,7)'>0]; % binary outcome (0: best, 1:worst)
        
        Y_PL = [Y_PL;y1;y3];
        U_PL = [U_PL;u1_PL;u3_PL];
        
        %--------------------------------------------------------------
        % FORMATING DATA FOR Q_LEARNING LEARNING
        %--------------------------------------------------------------
        % NOTES on how data should be formated
        % y: binary choice
        % u:
        %   - subject choice (binary)
        %   - reward for action 0
        %   - reward for action 1
        
        % outputs
        % experimenter data
          
        
        
          u1_QL = [y1;... % binary choice (0: correct, 1:incorrect)
     Outcomes(1, (data(i_cond1,7)'<0)+1 )];... % real outcome (0: best, 1:worst)
%     Outcomes(1, (data(i_cond1,7)'>0)+1 )];
 
           u3_QL = [y3;... % binary choice (0: correct, 1:incorrect)
     Outcomes(2, (data(i_cond3,7)'<0)+1 )];... % real outcome (0: best, 1:worst)
%     Outcomes(2, (data(i_cond3,7)'>0)+1 )]; 
 
        
        Y_QL = [Y_QL;y1;y3];
        U_QL = [U_QL;u1_QL;u3_QL];
          
        %=== second, let's consider the two conditions together : gains and loss.       
        
    end
    
    IsYout_QL = 0*Y_QL;
    IsYout_PL = 0*Y_PL;

    
    % Saving data after inversion
    clear data;
    data.Y_QL = Y_QL;
    data.Y_PL = Y_PL;
    data.U_QL = U_QL;
    data.U_PL = U_PL;
    data.IsYout_QL = IsYout_QL;
    data.IsYout_PL = IsYout_PL;
    data.i_subject = i_subject; % real index of subject (as in stefano's study)
    data.ind.cond1 = [1,3]; % which output indices correspond to condition 1
    data.ind.cond3 = [2,4]; % which output indices correspond to condition 3
    
    data_tache_mere{i_s}.data = data;
end





end

